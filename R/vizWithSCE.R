#' Visualize the integration of bulk DE results with scRNA-seq datasets
#'
#' This function assists with visualization of the integration
#' of bulk DE results tables with pre-processed scRNA-seq datasets available
#' on Bioconductor. After a user has picked a scRNA-seq dataset
#' using \code{integrateWithSingleCell} from the DESeq2 package,
#' they will have a results table, DESeqDataSet, and a SingleCellExperiment.
#' After annotating the SingleCellExperiment (see vignette), they can
#' provide that list of object to this function.
#' This function then produces a violin plot of the expression of
#' genes of interest from the bulk DE analysis across annotated cell types
#' in the single-cell dataset.
#'
#' @param dat the output of integrateWithSingleCell, a list with
#' results table (res), DESeqDataSet (dds), and SingleCellExperiment (sce)
#' @param which which gene to use for the violin plot.
#' This is represented as a number that corresponds to the lowest p-value
#' e.g. 1 for the lowest p-value, 2 for the second lowest p-value, etc.
#'
#' @return a violin plot (see Description)
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
