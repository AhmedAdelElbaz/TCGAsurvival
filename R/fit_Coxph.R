#' Fit a CoxPH Model on list of Genes
fit_Coxph <- function(Normalized_df , variables){

  d=rlang::parse_expr(paste0("survival::Surv(overall_survival, deceased)~",(paste(variables, collapse=" + "))))
  fit <- survival::coxph(as.formula(d) , data=Normalized_df)

return(fit)
}
