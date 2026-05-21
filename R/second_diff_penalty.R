second_diff_penalty <- function(k) {
  if (k <= 2) return(diag(1e-8, k))
  D2 <- diff(diag(k), differences = 2)
  t(D2) %*% D2 + diag(1e-8, k)
}
