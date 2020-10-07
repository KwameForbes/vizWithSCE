#' Visualize the integration of bulk DE results with Bioconductor single-cell RNA-seq datasets
#'
#' @description A package that assists with the visualization of the integration of bulk DE results tables with pre-processed scRNA-seq datasets available on Bioconductor, for downstream visualization tasks. After the user has pick a scRNA-seq dataset from integrateWithSingleCell.... The output of the function is...
#'
#'
#' @param dat this is a variable of integratewithSingleCell
#' @param which which gene to use for the violin plot
#'
#' @return nothing
#'
#' @examples
#'
#' \dontrun{
#' dat <- integrateWithSingleCell(res,dds)
#' vizWithSCE(dat)
#' }
#'
#' @author Michael Love
#' @author Kwame Forbes
#'
#' @export
vizWithSCE <- function(dat, which) {

  browser()
  
  stopifnot(all(names(dat) == c("res", "dds", "sce")))

  ## missing code for taking label and logcounts from sce

  # log counts? dat$sce => extracts log counts, but which gene?

  # gene according to 'which' => number of the gene in 'res' ranked by adjusted p-value

  o <- order(dat$res$padj)
  o[which] # top 'which' gene
  rownames(res)[o[which]] # name of the top 'which' gene

  dat$sce # this has log counts => index it by the name of the gene with the top 'which' adjusted p-value
  
  label <- colLabels(dat$sce)
  
  df <- data.frame(label, logcounts)
  
  ggplot(df, aes(label,logcounts)) +
    geom_violin(scale="width") +
    geom_sina(scale="width", alpha=.5)
}
