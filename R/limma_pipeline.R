#'Limma Function For Normalization
#'
#' @param object Summerizedexperimental object
#' @param condition_variable definition feature exported from colData from TCGA object
#' @param reference_group To decide which one of the values of your conditional variable is going to be the reference group
#' @param Norm.Method Normalization method , method=c("TMM","RLE","upperquartile","none")
#' @param plot a logical variable, if TRUE plot the mean variance scatter plot
#' @param rm.zeros a logical variable, if TRUE removes genes that have zero expression accross all samples before normalization
#' @return a list of voom: A complex object   contains the TMM+voom normalized data;
#'         eBayes: A complex object  contains the the fitted model plus a number of statistics related to each of the probes
#'         topTable: A simple table contains  the  differentially expressed genes sorted by p-value
#' @export
#' @examples TCGAnormalize_limma(TCGAsurvival::Brain_TCGA_GBM_Object , condition_variable= 'definition',reference_group="Solid Tissue Normal")
TCGAnormalize_limma = function(
    object,
    condition_variable= 'definition',
    reference_group="Solid Tissue Normal",
    Norm.Method = "TMM",
    plot =TRUE,
    rm.zeros = TRUE){
#Design the factor
design_factor =SummarizedExperiment::colData(object)[, condition_variable, drop=T]

group = factor(design_factor)
if(!is.null(reference_group)){group = relevel(group, ref=reference_group)}
#Matrix model for fitting
design = model.matrix(~ group)
#Create DGElist
dge =edgeR::DGEList(counts= SummarizedExperiment::assay(object),
                    samples=SummarizedExperiment::colData(object)[,],
                    genes=as.data.frame(SummarizedExperiment::rowData(object)),
                    remove.zeros = rm.zeros)

#filtering
keep =edgeR::filterByExpr(dge,design)
dge = dge[keep,,keep.lib.sizes=FALSE]
rm(keep)

# Normalization (TMM followed by voom)
dge =edgeR::calcNormFactors(dge, method = Norm.Method)#method=c("TMM","RLE","upperquartile","none")
v =limma::voom(dge, design, plot=plot)

# Fit model to data given design
fit =limma::lmFit(v, design)
fit = limma::eBayes(fit)

# Show top genes
topGenes =limma::topTable(fit, coef=ncol(design), number=dim(dge$genes)[1], sort.by="p")

return(
  list(
      voomObj=v, # normalized data
      fit=fit, # fitted model and statistics
      topGenes=topGenes # the differentially expressed genes ordered
    )
  )
}

