% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TCGAtumorSiteSelect.R
\name{TCGAtumorSiteSelect}
\alias{TCGAtumorSiteSelect}
\title{Prepare Object of Tumor Primary Site}
\usage{
TCGAtumorSiteSelect(
  tumor.site,
  sample.type = "Primary Tumor",
  project = "TCGA",
  DataCat = "Transcriptome Profiling",
  DataType = "Gene Expression Quantification",
  ExperStrategy = "RNA-Seq",
  nT = "all",
  nN = "all",
  S = 99,
  filepath = getwd()
)
}
\arguments{
\item{tumor.site}{primary tumor site.}

\item{sample}{type: by default "primary site".}

\item{file}{path: path for saving the summerizedExperiment object}
}
\value{
A summarizedExperiment
}
\description{
Prepare Object of Tumor Primary Site
}
\examples{
object <- TCGAtumorSiteSelect(tumor.site = "Kidney" , project = "TCGA", DataCat = "Transcriptome Profiling", DataType = "Gene Expression Quantification", ExperStrategy = "RNA-Seq", nT = 30 , nN = "all")
}
\author{
Ahmed Elbaz
}
