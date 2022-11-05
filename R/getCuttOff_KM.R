#' Get Cutoff With KM method
getCutoff_KM <- function(df){
  Q1 <- as.numeric(quantile(df$voomObjE)[2])
  Q3 <- as.numeric(quantile(df$voomObjE)[4])
  IQR <- as.vector(df$voomObjE[which(df$voomObjE >= Q1 & df$voomObjE < Q3)])

  y = c()
  v = c()
  df$deceased = df$vital_status == "Dead"
  df$overall_survival = as.numeric(ifelse(df$deceased,
                                          df$days_to_death,
                                          df$days_to_last_follow_up))

  for (V in IQR){
    df = df%>%mutate('gene_expr' = ifelse(df$`voomObjE` > V , "High" , "Low"))
    fit =survival::survfit(survival::Surv(overall_survival, deceased) ~ gene_expr , data = df)
    x <- survminer::surv_pvalue(fit,data = df)
    y = c(y , x$pval)
    v = c(v, V)
  }
  S = p.adjust(y)
  df <- data.frame(Value = v , P_value = y , adj.P_value = S)%>%dplyr::arrange(adj.P_value, P_value)


  return(df[1,1])
}

