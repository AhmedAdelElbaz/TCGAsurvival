#' Create a data frame for Survival Analysis
#'
#' @param object a SummerizedExperimental object.
#' @param limma_res Object by limma package
#' @param method the method to detect UP and DOWN expression level of certain gene according to soecific cut off, method = ("quantiles","median","KM_","CoxPH_")
#' @param columns.to.check character or vector of characters to be extracted from samples data.
#'
#' @return a data frame contains the gene expression row values, vital status, days to last follow up, normalized gene expression values and the level of expression as well as other columns specified.
#' @export
#' @examples GeneDataExtract(object = TCGAsurvival::Brain_TCGA_GBM_Object , Gene.list = c("FBP1" , "MAPK1") , limma_res = TCGAsurvival::Brain_limma_res , method = "KM_")
GeneDataExtract <- function(object,
                            Gene.list = Gene.list,
                            limma_res,
                            method = "KM_", #"quantiles"#"median" #CoxPH_"
                            plot = TRUE,
                            adj.P.value_Cutoff = 0.05,
                            logFC_Cutoff = 1,
                            columns.to.check = c()
)
{
  #Convert gene name into Ensebml ID
Gl <- c()
for (Gene in Gene.list){
  Gene = ifelse(stringr::str_detect(Gene,'ENSG'),
                Gene,
                ifelse(Gene %in% limma_res$voomObj$gene$gene_name,
                       get_row.name(limma_res$voomObj$genes,
                                    'gene_name',
                                    Gene),
                       stop(stringr::str_glue("ERROR: \n Selected Gene {Gene.list} may be Excluded due to low expression values")))
  )
  Gl = c(Gl, Gene)
#Density Plot for Gene expression across Tumor and Normal Samples
  if (plot){
    gExpM <- as.data.frame(SummarizedExperiment::assay(object))[Gene,]%>%tidyr::gather('sample','count')

    g <- ggplot2::ggplot(gExpM , ggplot2::aes(x = count))
    print(g +
            ggplot2::geom_density(fill = "#69b3a2", color = "#e9ecef", alpha=0.8)+
            ggplot2::ggtitle(stringr::str_glue('{Gene} Expr. Before Filteration'))+
            ggplot2::labs(y = "density", x = stringr::str_glue('{Gene}')))
  }
}
print(Gene.list)
gNames <- Gene.list
Gene.list = Gl
print(Gene.list)

  ##Step 3: Extract the selected genes for Survival analysis
  Normalized_df = extractGeneData(Query_Object = object ,
                                  Gene.list  = Gene.list,
                                  columns.to.check = c("vital_status" ,
                                                       "days_to_last_follow_up" ,
                                                       "definition",
                                                       "days_to_death",
                                                       "barcode", columns.to.check))



df_list = hash::hash()
  #Up and Down expressed genes
for (Gene in Gene.list){
  Normalized_df = Normalized_df%>%dplyr::mutate('voomObjE'
                                                = limma_res$voomObj$E[Gene,Normalized_df$barcode])

  #plot after Normalization
  g <- ggplot2::ggplot(Normalized_df, ggplot2::aes(x = voomObjE))
  print(g +
          ggplot2::geom_density(fill = "blue", color = "blue", alpha=0.5)+
          ggplot2::ggtitle(stringr::str_glue('{Gene} Expr. After Filteration'))+
          ggplot2::labs(y = "density", x = stringr::str_glue('{Gene}')))
  ####
  #Calculate the threshold for up and down
  ###

  #calculate with quantiles
  if (method == "quantiles"){
    if (limma_res$topGenes[Gene,"adj.P.Val"] < adj.P.value_Cutoff){
    #######################
      if(limma_res$topGenes[Gene,"logFC"] > logFC_Cutoff){
        Normalized_df = Normalized_df%>%mutate('gene_expr' = ifelse(Normalized_df$`voomObjE` >= quantile(Normalized_df$`voomObjE`)[2],ifelse(Normalized_df$`voomObjE` >= quantile(Normalized_df$`voomObjE`)[4], "Severly High", "Moderately High"),"Mild High"))

        print(stringr::str_glue('LogFC {Gene} = {limma_res$topGenes[Gene,"logFC"]}'))
        print(stringr::str_glue('adj.P.Val {Gene} = {limma_res$topGenes[Gene,"adj.P.Val"]}'))

      }else if (limma_res$topGenes[Gene,"logFC"] < -logFC_Cutoff){
        Normalized_df = Normalized_df%>%mutate('gene_expr' = ifelse(Normalized_df$`voomObjE` >= quantile(Normalized_df$`voomObjE`)[2],ifelse(Normalized_df$`voomObjE` >= quantile(Normalized_df$`voomObjE`)[4],"Mild Down" , "Moderately Down"), "Severly Down"))

        print(stringr::str_glue('LogFC {Gene} = {limma_res$topGenes[Gene,"logFC"]}'))
        print(stringr::str_glue('adj.P.Val {Gene} = {limma_res$topGenes[Gene,"adj.P.Val"]}'))
      }else{
        Normalized_df = Normalized_df%>%mutate('gene_expr' = ifelse(Normalized_df$`voomObjE` >= quantile(Normalized_df$`voomObjE`)[2],ifelse(Normalized_df$`voomObjE` >= quantile(Normalized_df$`voomObjE`)[4], "High", " Moderate"),"Low"))
        print(stringr::str_glue('LogFC {Gene} = {limma_res$topGenes[Gene,"logFC"]}'))
        print(stringr::str_glue('adj.P.Val {Gene} = {limma_res$topGenes[Gene,"adj.P.Val"]}'))
      }
  ########################
    }else{
      Normalized_df = Normalized_df%>%mutate('gene_expr' = ifelse(Normalized_df$`voomObjE` >= quantile(Normalized_df$`voomObjE`)[2],ifelse(Normalized_df$`voomObjE` >= quantile(Normalized_df$`voomObjE`)[4], "High", " Moderate"),"Low"))
      print(stringr::str_glue('LogFC {Gene} = {limma_res$topGenes[Gene,"logFC"]}'))
      print(stringr::str_glue('adj.P.Val {Gene} = {limma_res$topGenes[Gene,"adj.P.Val"]}'))
    }
    #Median method
  }else if (method == "median"){
    if (limma_res$topGenes[Gene,"adj.P.Val"] < adj.P.value_Cutoff){
    ############################
      if(limma_res$topGenes[Gene,"logFC"] > logFC_Cutoff){
        Normalized_df = Normalized_df%>%mutate('gene_expr' = ifelse(Normalized_df$`voomObjE` > getCutoff(Normalized_df), "Severly Up" , "Up"))
        print(stringr::str_glue('LogFC {Gene} = {limma_res$topGenes[Gene,"logFC"]}'))
        print(stringr::str_glue('adj.P.Val {Gene} = {limma_res$topGenes[Gene,"adj.P.Val"]}'))
      }else if (limma_res$topGenes[Gene,"logFC"] < -logFC_Cutoff){
        Normalized_df = Normalized_df%>%mutate('gene_expr' = ifelse(Normalized_df$`voomObjE` < getCutoff(Normalized_df), "Severly Down" , "Down"))
        print(stringr::str_glue('LogFC {Gene} = {limma_res$topGenes[Gene,"logFC"]}'))
        print(stringr::str_glue('adj.P.Val {Gene} = {limma_res$topGenes[Gene,"adj.P.Val"]}'))
      }else{
        Normalized_df = Normalized_df%>%mutate('gene_expr' = ifelse(Normalized_df$`voomObjE` > getCutoff(Normalized_df), "Higher than median" , "Lower than median"))
        print(stringr::str_glue('LogFC {Gene} = {limma_res$topGenes[Gene,"logFC"]}'))
        print(stringr::str_glue('adj.P.Val {Gene} = {limma_res$topGenes[Gene,"adj.P.Val"]}'))
      }
#############################
    }else{
      Normalized_df = Normalized_df%>%mutate('gene_expr' = ifelse(Normalized_df$`voomObjE` > getCutoff(Normalized_df), "Higher than median" , "Lower than median"))
      print(stringr::str_glue('LogFC {Gene} = {limma_res$topGenes[Gene,"logFC"]}'))
      print(stringr::str_glue('adj.P.Val {Gene} = {limma_res$topGenes[Gene,"adj.P.Val"]}'))

    }
  }else if (method == "KM_"){
    if (limma_res$topGenes[Gene,"adj.P.Val"] < adj.P.value_Cutoff){
      ############################
      if(limma_res$topGenes[Gene,"logFC"] > logFC_Cutoff){
        Normalized_df = Normalized_df%>%mutate('gene_expr' = ifelse(Normalized_df$`voomObjE` > getCutoff_KM(Normalized_df), "Severly Up" , "Up"))
        print(stringr::str_glue('LogFC {Gene} = {limma_res$topGenes[Gene,"logFC"]}'))
        print(stringr::str_glue('adj.P.Val {Gene} = {limma_res$topGenes[Gene,"adj.P.Val"]}'))
      }else if (limma_res$topGenes[Gene,"logFC"] < -logFC_Cutoff){
        Normalized_df = Normalized_df%>%mutate('gene_expr' = ifelse(Normalized_df$`voomObjE` < getCutoff_KM(Normalized_df), "Severly Down" , "Down"))
        print(stringr::str_glue('LogFC {Gene} = {limma_res$topGenes[Gene,"logFC"]}'))
        print(stringr::str_glue('adj.P.Val {Gene} = {limma_res$topGenes[Gene,"adj.P.Val"]}'))
      }else{
        Normalized_df = Normalized_df%>%mutate('gene_expr' = ifelse(Normalized_df$`voomObjE` > getCutoff_KM(Normalized_df), "High" , "Low"))
        print(stringr::str_glue('LogFC {Gene} = {limma_res$topGenes[Gene,"logFC"]}'))
        print(stringr::str_glue('adj.P.Val {Gene} = {limma_res$topGenes[Gene,"adj.P.Val"]}'))
      }
      #############################
    }else{
      Normalized_df = Normalized_df%>%mutate('gene_expr' = ifelse(Normalized_df$`voomObjE` > getCutoff_KM(Normalized_df), "High" , "Low"))
      print(stringr::str_glue('LogFC {Gene} = {limma_res$topGenes[Gene,"logFC"]}'))
      print(stringr::str_glue('adj.P.Val {Gene} = {limma_res$topGenes[Gene,"adj.P.Val"]}'))

    }
  }else if (method == "CoxPH_"){
    if (limma_res$topGenes[Gene,"adj.P.Val"] < adj.P.value_Cutoff){
      ############################
      if(limma_res$topGenes[Gene,"logFC"] > logFC_Cutoff){
        Normalized_df = Normalized_df%>%mutate('gene_expr' = ifelse(Normalized_df$`voomObjE` > getCutoff_CoxPH(Normalized_df), "Severly Up" , "Up"))
        print(stringr::str_glue('LogFC {Gene} = {limma_res$topGenes[Gene,"logFC"]}'))
        print(stringr::str_glue('adj.P.Val {Gene} = {limma_res$topGenes[Gene,"adj.P.Val"]}'))
      }else if (limma_res$topGenes[Gene,"logFC"] < -logFC_Cutoff){
        Normalized_df = Normalized_df%>%mutate('gene_expr' = ifelse(Normalized_df$`voomObjE` < getCutoff_CoxPH(Normalized_df), "Severly Down" , "Down"))
        print(stringr::str_glue('LogFC {Gene} = {limma_res$topGenes[Gene,"logFC"]}'))
        print(stringr::str_glue('adj.P.Val {Gene} = {limma_res$topGenes[Gene,"adj.P.Val"]}'))
      }else{
        Normalized_df = Normalized_df%>%mutate('gene_expr' = ifelse(Normalized_df$`voomObjE` > getCutoff_CoxPH(Normalized_df), "High" , "Low"))
        print(stringr::str_glue('LogFC {Gene} = {limma_res$topGenes[Gene,"logFC"]}'))
        print(stringr::str_glue('adj.P.Val {Gene} = {limma_res$topGenes[Gene,"adj.P.Val"]}'))
      }
      #############################
    }else{
      Normalized_df = Normalized_df%>%mutate('gene_expr' = ifelse(Normalized_df$`voomObjE` > getCutoff_CoxPH(Normalized_df), "High" , "Low"))
      print(stringr::str_glue('LogFC {Gene} = {limma_res$topGenes[Gene,"logFC"]}'))
      print(stringr::str_glue('adj.P.Val {Gene} = {limma_res$topGenes[Gene,"adj.P.Val"]}'))

    }
  }


df_list[[Gene]] = Normalized_df
}
#Change Column names
if (length(Gene.list) > 1){

  Normalized_df = df_list[[Gene.list[1]]]

  colnames(Normalized_df)[which(colnames(Normalized_df) == 'voomObjE')] = stringr::str_glue('{Gene.list[1]}_norm.expr.')

  colnames(Normalized_df)[which(colnames(Normalized_df) == 'gene_expr')] = stringr::str_glue('{Gene.list[1]}_Level')

  for (Gene in Gene.list[2:length(Gene.list)]){
    Normalized_df = cbind(Normalized_df ,
                          df_list[[Gene]][['voomObjE']],
                          df_list[[Gene]][['gene_expr']])
    colnames(Normalized_df)[which(colnames(Normalized_df) == 'df_list[[Gene]][["voomObjE"]]')] = stringr::str_glue('{Gene}_norm.expr.')
    colnames(Normalized_df)[which(colnames(Normalized_df) == 'df_list[[Gene]][["gene_expr"]]')] = stringr::str_glue('{Gene}_Level')

  }
}else{
  Normalized_df = df_list[[`Gene.list`]]
}

  return(Normalized_df)
}
