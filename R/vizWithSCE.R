#' Visualize the integration of bulk DE results with Bioconductor single-cell RNA-seq datasets
#'
#' @description A package that assists with the visualization of the integration of bulk DE results tables with pre-processed scRNA-seq datasets available on Bioconductor, for downstream visualization tasks. After the user has pick a scRNA-seq dataset from integrateWithSingleCell.... The output of the function is...
#'
#'
#' @param dat this is a variable of integratewithSingleCell
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
vizWithSCE <- function(dat) {
  stopifnot(all(names(dat) == c("res", "dds", "sce")))
  ggplot(dat, aes(label,logcounts)) +
    geom_violin(scale="width") +
    geom_sina(scale="width", alpha=.5) +
    ggtitle("Cd52 expression across clusters")

  plotColData(sce, y=I(logcts), x="label")
}
integrateWithSingleCell(res,dds)
getwd()
