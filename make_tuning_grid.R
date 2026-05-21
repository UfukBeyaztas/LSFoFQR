make_tuning_grid <- function(n_train,
                             base_tuning = default_tuning(n_train),
                             smooth_mult = c(0.5, 1.0, 2.0),
                             group_mult  = c(0.5, 1.0, 2.0),
                             active_tol_grid = NULL,
                             eta_grid = NULL,
                             active_method = "fixed_or_gap",
                             min_gap = 0.02) {
  if (is.null(active_tol_grid)) {
    active_tol_grid <- sort(unique(pmin(0.55, pmax(0.03,
                                                   base_tuning$active_tol + c(-0.08, 0, 0.08)
    ))))
  }
  if (is.null(eta_grid)) eta_grid <- base_tuning$eta_adapt
  
  raw_grid <- expand.grid(
    smooth_mult = smooth_mult,
    group_mult  = group_mult,
    active_tol  = active_tol_grid,
    eta_adapt   = eta_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  out <- vector("list", nrow(raw_grid))
  for (rr in seq_len(nrow(raw_grid))) {
    sm <- raw_grid$smooth_mult[rr]
    gm <- raw_grid$group_mult[rr]
    out[[rr]] <- list(
      lam_init_b = base_tuning$lam_init_b * sm,
      lam_init_a = base_tuning$lam_init_a * sm,
      lam_sp_b   = base_tuning$lam_sp_b   * sm,
      lam_sp_a   = base_tuning$lam_sp_a   * sm,
      lam_sp_g   = base_tuning$lam_sp_g   * gm,
      lam_rf_b   = base_tuning$lam_rf_b   * sm,
      lam_rf_a   = base_tuning$lam_rf_a   * sm,
      active_tol = raw_grid$active_tol[rr],
      eta_adapt  = raw_grid$eta_adapt[rr],
      active_method = active_method,
      min_gap = min_gap
    )
  }
  out
}