ssolve <- function(A, b, ridge = 1e-8) {
  A <- as.matrix(A)
  p <- ncol(A)
  out <- tryCatch(
    as.numeric(solve(A + diag(ridge, p), b)),
    error = function(e) as.numeric(ginv(A + diag(ridge, p)) %*% b)
  )
  out
}
