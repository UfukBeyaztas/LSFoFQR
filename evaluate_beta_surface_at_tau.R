evaluate_beta_surface_at_tau <- function(beta_vec, bases, tau_index) {
  Kt <- ncol(bases$Bt)
  Ks <- ncol(bases$Bs)
  Ku <- ncol(bases$Bu)
  
  beta_arr <- array(beta_vec, c(Kt, Ks, Ku))
  out <- matrix(0.0, nrow(bases$Bt), nrow(bases$Bs))
  
  for (k in seq_len(Ku)) {
    out <- out + as.matrix(bases$Bt %*% beta_arr[, , k] %*% t(bases$Bs)) * bases$Bu[tau_index, k]
  }
  
  out
}