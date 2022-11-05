#' Get Row Name
get_row.name = function(df,getter,x){
  df = df%>%tidyr::as_tibble()%>%mutate("rn" = row.names(df))%>%select(c(rn,getter))%>%filter(`gene_name`==x)
  y = df[,"rn"]
  return(as.character(y))
}
