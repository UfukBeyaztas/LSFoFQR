make_ridge_warm_start <- function(design, Q_b, Q_a, ridge = 1e-8) {
  A <- cbind(design$X_beta, design$X_alpha)
  Q <- bdiag(Q_b, Q_a)
  H <- crossprod(A) / design$N + Q
  rhs <- as.numeric(crossprod(A, design$y) / design$N)
  ssolve(H, rhs, ridge = ridge)
}
