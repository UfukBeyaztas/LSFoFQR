predict.sfqfofr <- function(object, Xnew, ...) {
  
  fit   <- object$fit_r
  grids <- object$grids
  bases <- object$bases
  integration <- object$integration
  
  R_new <- compute_score_matrix(
    Xmat = Xnew,
    Bt = bases$Bt,
    t_grid = grids$t,
    integration = integration
  )
  
  Kt <- ncol(bases$Bt)
  Ks <- ncol(bases$Bs)
  Ku <- ncol(bases$Bu)
  
  ns <- length(grids$s)
  nu <- length(grids$tau)
  n  <- nrow(Xnew)
  
  beta_arr  <- array(fit$beta,  c(Kt, Ks, Ku))
  alpha_arr <- array(fit$alpha, c(ncol(bases$Bsa), ncol(bases$Bua)))
  
  out <- array(0.0, c(n, ns, nu))
  
  for (r in seq_len(nu)) {
    Theta_r <- matrix(0.0, Kt, Ks)
    
    for (k in seq_len(Ku)) {
      Theta_r <- Theta_r + beta_arr[, , k] * bases$Bu[r, k]
    }
    
    alpha_r <- as.numeric(bases$Bsa %*% alpha_arr %*% bases$Bua[r, ])
    
    out[, , r] <- sweep(
      R_new %*% Theta_r %*% t(bases$Bs),
      2L,
      alpha_r,
      "+"
    )
  }
  
  out
}
