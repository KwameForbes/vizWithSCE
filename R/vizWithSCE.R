install.packages("devtools")
library(devtools)
devtools::install_github
install_github("mikelove/DESeq2")
library(DESeq2)
packageVersion("DESeq2")


plotter <- function(dat) {
  stopifnot(all(names(dat) == c("res", "dds", "ans")))
  plot(dat$res, dat$dds)
}
