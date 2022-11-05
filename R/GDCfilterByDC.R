#' GDC filter based on data category
#'
#' @param data.category : The data category entered by the user.
#' @return list of each project that has the desired data category entered by the user.
#' @details It enters to GDC project , retrieve the size of each project. Then, it shows the data category of each project. Finally, it  filters the data based on the desired data category.
GDCfilterByDC <- function(data.category){

  GDCprojects = TCGAsurvival::GDC_projects

  GDCprojectsize <- function(){

    project.sizes = list()
    for (i in 1:length(GDCprojects)){
      project.sizes[[names(GDCprojects)[i]]] = GDCprojects[[names(GDCprojects)[i]]]["file_size"]
    }
    return(project.sizes)
  }


  projects_sizes.df = as.data.frame(GDCprojectsize())
  project_size.data_cat = data.frame(id = names(GDCprojectsize()) , size = t(projects_sizes.df) )

  project.data.category = list()
  for (i in 1:length(names(GDCprojects))){
    project.data.category[[names(GDCprojects)[i]]] = as.data.frame(GDCprojects[[names(GDCprojects)[i]]][2])$data_categories.data_category
  }


  filtered_projects_by.data.category = list()
  for (i in 1:length(names(project.data.category))){
    for (j in 1:length(project.data.category[[i]])){
      if (project.data.category[[i]][j] == data.category){
        filtered_projects_by.data.category[names(project.data.category)[i]] = project.data.category[[i]][j]
      }else{
        next
      }
    }
  }
  project_size.data_cat["found data_category"] = "not found"


  for (j in 1:length(names(filtered_projects_by.data.category))){
    for(i in 1:nrow(project_size.data_cat)){
      if (names(filtered_projects_by.data.category)[j] == project_size.data_cat$id[i]){
        project_size.data_cat$`found data_category`[i] = data.category
      }else{
        next
      }
    }
  }
  row.names(project_size.data_cat) = project_size.data_cat$`id`
  print(table(project_size.data_cat$`found data_category`))
  return(list(filteredGDC = project_size.data_cat,
              AllGDC = project.data.category,
              GDCprojSummary = GDCprojects))
}
