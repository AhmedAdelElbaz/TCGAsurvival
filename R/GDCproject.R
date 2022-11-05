#' Function retrieves All GDC project summary
#'
#'@details GDC project function return list of project IDs including summary of each project. This summary includes data category, case count and file count.
#' @return environment with the summary of each project in the GDC data base.
#' @export
#' @examples GDCproject()
GDCproject <- function(){
  GDC_projects = TCGAbiolinks::getGDCprojects()
  GDCprojects = hash::hash()
  for (id in GDC_projects$id){
    GDCprojects[[id]] = TCGAbiolinks::getProjectSummary(id)
  }
  rm(id)
  return(GDCprojects)
}
