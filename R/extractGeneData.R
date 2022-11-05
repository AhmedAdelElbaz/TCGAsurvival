#' Extract Gene Data
extractGeneData <- function(Query_Object,
                            Gene.list = genes$gene_id ,
                            sample.type = 'Primary solid Tumor',
                            columns.to.check = names(clinical)){
  #Select Tumor
  clinical <- Query_Object@colData%>%tidyr::as_tibble()%>%dplyr::filter(definition == sample.type)
  gExpM <- SummarizedExperiment::assay(Query_Object)%>%tidyr::as_tibble()%>%select(clinical$barcode)
  expression_df <- as.data.frame(gExpM)

  expression_df <-dplyr::mutate(expression_df ,
                                'Ensemble ID' = row.names(assay(Query_Object)))
  tidy_gene_expression <-tidyr::as_tibble(expression_df)%>%tidyr::gather(sample, expression , -'Ensemble ID')%>%tidyr::spread(`Ensemble ID`,expression)%>%select(sample , Gene.list)
  subsetted_clinical =tidyr::as_tibble(clinical)%>%dplyr::select('barcode' , columns.to.check)
  genes <- SummarizedExperiment::rowData(Query_Object)

  total <- integrate_files(to_dataframe =  tidy_gene_expression ,
                           from_dataframe =  subsetted_clinical ,
                           key_to =  'sample' ,
                           key_from =  'barcode' ,
                           columns.to.check =  columns.to.check)

  return(total)
}
