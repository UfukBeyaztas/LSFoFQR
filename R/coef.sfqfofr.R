coef.sfqfofr <- function(object, tau_index = NULL, ...) {
  
  fit_r <- object$fit_r
  bases <- object$bases
  nu <- length(object$grids$tau)
  
  nt <- nrow(bases$Bt)
  ns <- nrow(bases$Bs)
  
  arr <- array(0.0, c(nt, ns, nu))
  
  for (r in seq_len(nu)) {
    arr[, , r] <- evaluate_beta_surface_at_tau(fit_r$beta, bases, r)
  }
  
  if (!is.null(tau_index)) {
    return(arr[, , tau_index, drop = FALSE])
  }
  
  arr
}
