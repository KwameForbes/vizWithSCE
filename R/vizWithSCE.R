#' Visualize the integration of bulk DE results with scRNA-seq datasets
#'
#' @description A package that assists with the visualization of the integration
#' of bulk DE results tables with pre-processed scRNA-seq datasets available
#' on Bioconductor, for downstream visualization tasks. After the user has pick
#' a scRNA-seq dataset from integrateWithSingleCell....
#' The output of the function is...
#'
#'
#' @param dat the output of integrateWithSingleCell, a list with
#' results table (res), DESeqDataSet (dds), and SingleCellExperiment (sce)
#' @param which which gene to use for the violin plot.
#' This is represented as a number that corresponds to the lowest p-value
#' e.g. 1 for the lowest p-value, 2 for the second lowest p-value, etc.
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
#' @author Kwame Forbes
#' @author Michael Love 
#'
#' @importFrom ggforce geom_sina
#' 
#' @export
vizWithSCE <- function(dat, which) {

  #stopifnot(all(names(dat) == c("res", "dds", "sce")))

  ## missing code for taking label and logcounts from sce
  # log counts? dat$sce => extracts log counts, but which gene?
  # gene according to 'which' => number of the gene in 'res' ranked by p-value

  res <- dat$res
  dds <- dat$dds
  sce <- dat$sce

  o <- order(res$pvalue)
  o[which] # top 'which' gene
  gene <-rownames(res)[o[which]] # name of the top 'which' gene

  stopifnot(gene %in% rownames(sce))

  # this has log counts => index it by the name of the gene
  # with the top 'which' adjusted p-value
  sce 

  label <- colLabels(sce)

  df <- data.frame(label=colLabels(sce), logcounts=logcounts(sce)[gene,])

  ggplot(df, aes(label,logcounts)) + geom_violin(scale="width") +
    ggforce::geom_sina(scale="width", alpha=.5)
}
