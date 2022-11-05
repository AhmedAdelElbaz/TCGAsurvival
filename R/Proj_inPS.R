#' Each Primary Tumor Site & The Corresponding Projects
#'
#' @return: Dict. or Environment, each key represents a primary tumor site and the values are the coressponding projects
#' @export
#' @examples PS <- Projects_in_Primary_sites()
Projects_in_Primary_sites <- function(){
  PS_dict = hash::hash()
  GDC_projects =TCGAbiolinks::getGDCprojects()
  #Primary Sites
  extract <- function(list){
    X = c()
    for (i in 1:length(list)){
      x = list[[i]]
      X = append(X , x)
    }
    result = as.data.frame(table(X))
    names(result) = c("primary_sites", "count")
    return(result)
  }

  Primary_sites_in_projects = extract(GDC_projects$primary_site)
  for (ps in Primary_sites_in_projects$primary_sites){
    empty_list  = c()
    for (i in 1:length(GDC_projects$primary_site)){
      if (length(GDC_projects$primary_site[[i]]) != 0){

        for (j in 1:length(GDC_projects$primary_site[[i]])){

          if (ps == GDC_projects$primary_site[[i]][j]){
            empty_list = append(empty_list , GDC_projects$id[[i]])

          }
        }

      }

    }
    PS_dict[[ps]] = empty_list
  }

  return(PS_dict)
}

