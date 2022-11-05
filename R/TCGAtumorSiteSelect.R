#' Prepare Object of Tumor Primary Site
#'
#' @param tumor.site primary tumor site.
#' @param sample type: by default "primary site".
#' @param file path: path for saving the summerizedExperiment object
#' @return A summarizedExperiment
#' @export
#' @author Ahmed Elbaz
#' @examples object <- TCGAtumorSiteSelect(tumor.site = "Kidney" , project = "TCGA", DataCat = "Transcriptome Profiling", DataType = "Gene Expression Quantification", ExperStrategy = "RNA-Seq", nT = 30 , nN = "all")
TCGAtumorSiteSelect <- function(tumor.site ,
                                sample.type = 'Primary Tumor',
                                project = "TCGA",
                                DataCat = "Transcriptome Profiling",
                                DataType = "Gene Expression Quantification",
                                ExperStrategy = "RNA-Seq",
                                nT= "all",
                                nN = "all",
                                S = 99,
                                filepath = getwd()){

#Projects for Primary Site
psPid = TCGAsurvival::Primary_sites
#Select the Project
  ifelse(length(psPid[[tumor.site]][which(stringr::str_detect(psPid[[tumor.site]], project))]) == 0,
         stop("Selected Primary Site not included in TCGA projects, please choose another Primary site or use TCGAbiolinks to create query object"),
         print(tumor.site))

  for (PS in names(psPid)){
    if (tumor.site == PS){

      project_id = psPid[[tumor.site]][which(stringr::str_detect(psPid[[tumor.site]], project))]
      #Get the highest project with case counts
      n <- TCGAprojectCaseCount()[project_id,]%>%dplyr::arrange(dplyr::desc(count))
      project_id = row.names(n)[1]

#Create the Query
query = TCGAbiolinks::GDCquery(project = project_id,
                               data.category = DataCat,
                               data.type = DataType,
                               access = "open",
                               data.format = "TSV",
                               experimental.strategy = ExperStrategy,
      )

#filtering Primary Tumor and Normal samples only
nT = ifelse(nT == "all", dim(query$results[[1]]%>%filter(sample_type == sample.type))[1], nT)
nN = ifelse(nN == "all", dim(query$results[[1]]%>%filter(sample_type == "Solid Tissue Normal"))[1], nN)
#Set the required number of samples
barcode <- SampleSeclect(query,
                  sample.type = sample.type,
                  nT = nT,
                  nN = nN,
                  S = S)
#Reselect the samples
query = TCGAbiolinks::GDCquery(project = project_id,
                               data.category = DataCat,
                               data.type = DataType,
                               access = "open",
                               data.format = "TSV",
                               experimental.strategy = ExperStrategy,
                               barcode = barcode
)
  }#End of if condition
}#End of for loop
#Download samples
TCGAbiolinks::GDCdownload(query)
#Object
Query_Object = TCGAbiolinks::GDCprepare(query)
#Object Saving
fileName = stringr::str_glue('{filepath}/{tumor.site}_{project_id}_Object.RDS')
saveRDS(object = Query_Object , file = fileName , compress = FALSE)

return(Query_Object)
}
