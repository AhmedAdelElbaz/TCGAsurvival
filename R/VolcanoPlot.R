#' Plotting a Volcano PLot Display DEGs
#'
#' @param limma_res an object (output of TCGAnormalize_limma()). the object contains the normalized genes' expression values and the fitted model, as well as the genes' names.
#' @param logFC_Cutoff Selected cut off of the log Fold-Change
#' @param adj.p.Value_Cutoff adjusted P value cutoff to evaluate the differential gene expression significance.
#'
#' @return displays a Volcano plot with Up and Down differentially expressed Genes
#' @export
#' @examples volcPlot(limma_res = TCGAsurvival::Brain_limma_res, logFC_Cutoff = 2,adj.p.Value_Cutoff = 0.05)
volcPlot <- function(limma_res,
                     logFC_Cutoff = 2,
                     adj.p.Value_Cutoff = 0.05
                     #colorList = c("pink" , "black" , "blue")
                     ){
#Create a data frame for plotting
df <- data.frame(adj.pval =limma_res$topGenes$adj.P.Val,
                 logFC = limma_res$topGenes$logFC,
                 gene_name = limma_res$topGenes$gene_name)

df <- dplyr::mutate(df , 'gene_expr' = ifelse(df$adj.pval <= 0.05,
                                              ifelse(df$logFC > logFC_Cutoff , "UP" , ifelse(df$logFC < -logFC_Cutoff , "DOWN", "NO")), "NO"))
#Add Gene Names
df$gene_name <- ifelse(df$gene_expr == "NO" , NA , df$gene_name)
#Create ggplot2 object
g <- ggplot2::ggplot(data = df , ggplot2::aes(x = logFC , y = -log10(adj.pval), color = gene_expr , label = gene_name))
#Plotting
plot = g + ggplot2::geom_point() + ggrepel::geom_text_repel() +
  ggplot2::theme_minimal() +
  #scale_color_manual(values=colorList) +
  ggplot2::geom_vline(xintercept=c(-logFC_Cutoff, logFC_Cutoff), col="red") +
  ggplot2::geom_hline(yintercept=-log10(adj.p.Value_Cutoff), col="red")
  print(plot)

}







