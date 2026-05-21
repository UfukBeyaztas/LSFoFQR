coefficient_multiplicity <- function(groups, p) {
  mult <- tabulate(unlist(groups, use.names = FALSE), nbins = p)
  mult[mult == 0] <- 1
  mult
}
