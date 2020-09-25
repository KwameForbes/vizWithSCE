plotter <- function(dat) {
  stopifnot(all(names(dat) == c("res", "dds", "ans")))
  plot(dat$res, dat$dds)
}
