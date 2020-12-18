#' Visualize the integration of bulk DE results with Bioconductor single-cell RNA-seq datasets
#'
#' @description A package that assists with the visualization of the integration of bulk DE results tables with pre-processed scRNA-seq datasets available on Bioconductor, for downstream visualization tasks. After the user has pick a scRNA-seq dataset from integrateWithSingleCell.... The output of the function is...
#'
#'
#' @param dat this is a variable of integrateWithSingleCell
#' @param which which gene to use for the violin plot. This is represented as a number that corresponds to the lowest adjusted p-value
#' e.g. 1 for the lowest adjusted p-value, 2 for the second lowest adjusted p-value, 3 for the third, etc.
#'
#' @return nothing
#'
#' @examples
#'
#' \dontrun{
#' dat <- integrateWithSingleCell(res,dds)
#' vizWithSCE(dat, which=1)
#' }
#'
#' @author Michael Love
#' @author Kwame Forbes
#'
#' @export
vizWithSCE <- function(dat, which) {



  #stopifnot(all(names(dat) == c("res", "dds", "sce")))

  ## missing code for taking label and logcounts from sce
  # log counts? dat$sce => extracts log counts, but which gene?
  # gene according to 'which' => number of the gene in 'res' ranked by adjusted p-value

  res <- dat$res
  dds <- dat$dds
  sce <- dat$sce

  o <- order(res$padj)
  o[which] # top 'which' gene
  gene <-rownames(res)[o[which]] # name of the top 'which' gene

  stopifnot(gene %in% rownames(sce))

  sce # this has log counts => index it by the name of the gene with the top 'which' adjusted p-value

  label <- colLabels(sce)

  df <- data.frame(label=colLabels(sce), logcounts=logcounts(sce)[gene,])

  ggplot(df, aes(label,logcounts)) + geom_violin(scale="width")  +
  ggforce::geom_sina(scale="width", alpha=.5)
}
# vizWithSCE(dat, which=3)


