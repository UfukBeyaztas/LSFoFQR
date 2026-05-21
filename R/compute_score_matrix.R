compute_score_matrix <- function(Xmat,
                                 Bt,
                                 t_grid,
                                 integration = c("trapezoid", "simpson", "auto")) {
  
  integration <- match.arg(integration)
  
  Xmat <- as.matrix(Xmat)
  Bt   <- as.matrix(Bt)
  
  if (ncol(Xmat) != length(t_grid)) {
    stop("ncol(Xmat) must equal length(t_grid).")
  }
  
  if (nrow(Bt) != length(t_grid)) {
    stop("nrow(Bt) must equal length(t_grid).")
  }
  
  wt <- integration_weights(t_grid, method = integration)
  
  Xmat %*% (wt * Bt)
}
