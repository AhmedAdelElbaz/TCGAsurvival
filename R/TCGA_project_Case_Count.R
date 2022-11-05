#' TCGA_Project_Case_Count
TCGAprojectCaseCount <- function(){
  GDCprojects = TCGAsurvival::GDC_projects
  projCC <- data.frame(project = names(GDCprojects),
                       count = rep(0, length(GDCprojects)))
  for (p in names(GDCprojects)){
    x <- as.integer((GDCprojects[[p]]$data_categories%>%filter(data_category == 'Transcriptome Profiling'))$case_count)
    if (length(x) != 0){projCC[projCC$project == p,]$count = x}
  }
  row.names(projCC) = names(GDCprojects)
  return(projCC)
}
