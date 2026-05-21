tuning_list_to_row <- function(tuning, id = NA_integer_) {
  data.frame(
    id = id,
    lam_init_b = tuning$lam_init_b,
    lam_init_a = tuning$lam_init_a,
    lam_sp_b   = tuning$lam_sp_b,
    lam_sp_a   = tuning$lam_sp_a,
    lam_sp_g   = tuning$lam_sp_g,
    lam_rf_b   = tuning$lam_rf_b,
    lam_rf_a   = tuning$lam_rf_a,
    active_tol = tuning$active_tol,
    eta_adapt  = tuning$eta_adapt,
    active_method = tuning$active_method,
    min_gap = tuning$min_gap,
    stringsAsFactors = FALSE
  )
}