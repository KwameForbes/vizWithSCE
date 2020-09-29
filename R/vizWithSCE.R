install_github("mikelove/DESeq2")
library(DESeq2)
package Version("DESeq2)

plotter <- function(dat) {
  stopifnot(all(names(dat) == c("res", "dds", "ans")))
  plot(dat$res, dat$dds)
}
