#' Sample Select
#'
#' @param query: an object (output of TCGAtumorSiteSelect(...))
#' @param sample.type: Type of tumor
#' @param nT: number of Tumor samples
#' @param nN: number of Normal samples
#' @param S: set.seed function argument
#' @return list of samples to be downloaded from GDC database
SampleSeclect <- function(query,
                          sample.type = 'Primary Tumor',
                          nT ,
                          nN ,
                          S = S){
  set.seed(S)
  res <- query$results[[1]]
  dfTumor <- res%>%dplyr::filter(sample_type == sample.type)%>%sample_n(nT)

  dfNorm <- res%>%dplyr::filter(sample_type == 'Solid Tissue Normal')%>%sample_n(nN)
  barcode <- c(dfTumor$cases, dfNorm$cases)
  return(barcode)
}
