#'  Get Cutoff With Median method
getCutoff <- function(df){
  Q1 <- as.numeric(quantile(df$voomObjE)[2])
  Q3 <- as.numeric(quantile(df$voomObjE)[4])
  IQR <- as.vector(df$voomObjE[which(df$voomObjE >= Q1 & df$voomObjE < Q3)])

  y = c()
  v = c()

  for (V in IQR){
    Cf <- V
    d1 <- df%>%filter(voomObjE>Cf)%>%select(voomObjE)
    d2 <- df%>%filter(voomObjE<Cf)%>%select(voomObjE)
    s <- (wilcox.test(d1$voomObjE,
                      d2$voomObjE))
    y = c(y , s$p.value)
    v = c(v, V)
  }
  df <- data.frame(Value = v , P_value = y)%>%dplyr::arrange(P_value)


  return(df[1,1])
}
