#' PCA Plotting Funcion
#'
#' @param limma_res: an object, (output of TCGAnormalize_limma(...))
#' @param condition_variable: "definition" column in the limma_res object by default.
#' @param Pos: Position of the legend in the plot (default #topleft) #buttomleft #topright ..
#' @return: a data frame with all principal components.
#' @export
#'
#' @examples PCAplot(limma_res = TCGAsurvival::Brain_limma_res)
#'
PCAplot <- function(limma_res, condition_variable = "definition", Pos = "topleft"){
  group = factor(limma_res$voomObj$targets[, condition_variable])
  pca = prcomp(t(limma_res$voomObj$E))
  # Take PC1 and PC2 for the plot
  print(plot(pca$x[,1:2],col=group, pch=19))
  # include a legend for points
  legend(Pos, inset=.01, levels(group), pch=19, col=1:length(levels(group)))

  return(pca)
}

#This plot shows that the two sample groups (tumor tissue, and healthy tissue) have a well-separated RNA expression profile.


