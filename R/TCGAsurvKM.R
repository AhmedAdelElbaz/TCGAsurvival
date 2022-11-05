#' Fit a Survival Model For a List of Genes
#'
#' @param Normalized_df: a data frame contains the normalized gene (or list of genes) expression values, and its corresponding level of expression.
#' @param plot: logical. if TRUE the KM plot will be displayed.
#' @return For method "KM" returns fitted model including P Value of the model.
#' @export
TCGAsurv <- function(Normalized_df,
                       limma_res,
                       Gene.list,
                       method = "KM",#CoxPH
                       plot = TRUE)
{

  ##Step 4: Survival Analysis
  Normalized_df$deceased = Normalized_df$vital_status == "Dead"
  Normalized_df$overall_survival = as.numeric(ifelse(Normalized_df$deceased,
                                                     Normalized_df$days_to_death,
                                                     Normalized_df$days_to_last_follow_up))


  variables <- c()
  fit.list <- hash::hash()
  #gEnsm <- c()
  #logFC_list <- c()
  #adj.P.Value.logFC <- c()

#Convert from Ensemble to Gene names
for (Gene in Gene.list){
    Gene = ifelse(stringr::str_detect(Gene,'ENSG'),
                  Gene,
                  ifelse(Gene %in% limma_res$voomObj$gene$gene_name,
                         get_row.name(limma_res$voomObj$genes,
                                      'gene_name',
                                      Gene),
                         stop(stringr::str_glue("ERROR: \n Selected Gene {Gene.list} may be Excluded due to low expression values")))
    )

  #gEnsm <- c(gEnsm , Gene)
    colnames(Normalized_df)[which(colnames(Normalized_df) == stringr::str_glue('{Gene}_Level'))] = 'gene_expr'


#Fit the Model
if (method == "KM"){


  fit =survival::survfit(survival::Surv(overall_survival, deceased) ~ gene_expr , data = Normalized_df)

#print(stringr::str_glue("{Gene} survfit P.value = {survminer::surv_pvalue(fit)}"))

  if (plot){ploty <-survminer::ggsurvplot(fit, data=Normalized_df,
                                          pval=T,
                                          risk.table=T,
                                          risk.table.col="strata",
                                          title = stringr::str_glue('{Gene}_Level'),
                                          legend.title = "no. of Samples",
                                          xlab = "Time by Days"
                                          )
  print(ploty)

  }

  colnames(Normalized_df)[which(colnames(Normalized_df) == 'gene_expr')] = stringr::str_glue('{Gene}_Level')

  fit.list[[Gene]] = fit


}else if (method == "CoxPH"){

if (length(Gene.list) == 1){
  if ('voomObjE' %in% names(Normalized_df)){
    variables = 'voomObjE'
  }else{
    variables = c(variables , names(Normalized_df[which(stringr::str_detect(names(Normalized_df),pattern = stringr::str_glue("{Gene}_norm.expr.")))]))
  }

}else{
  variables = c(variables , names(Normalized_df[which(stringr::str_detect(names(Normalized_df),pattern = stringr::str_glue("{Gene}_norm.expr.")))]))

}

fit <- fit_Coxph(Normalized_df,
                 variables)


fit.list = summary(fit)
fit.list <- data.frame(Genes = Gene.list ,
                       #"Ensmble ID" = gEnsm,
                       #logFC_list,
                       #adj.pVal = adj.P.Value.logFC,
                       coef = fit.list$coefficients[,1][1:length(Gene.list)],
                       "exp(coef)" = fit.list$coefficients[,2][1:length(Gene.list)],
                       "se(coef)" = fit.list$coefficients[,3][1:length(Gene.list)],
                       z = fit.list$coefficients[,4][1:length(Gene.list)],
                       P.value = fit.list$coefficients[,5][1:length(Gene.list)])
row.names(fit.list) = seq(1,length(Gene.list))

}

}
  return(fit.list)
}


