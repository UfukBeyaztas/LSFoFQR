sfqfofr_tuning <- function(design_train, bases, ginfo, M_grp, tuning) {
  fit_p <- fit_pilot(
    design_train, bases, ginfo,
    lambda_beta = tuning$lam_init_b,
    lambda_alpha = tuning$lam_init_a,
    ridge_beta = 2e-5,
    ridge_alpha = 1e-5,
    rough_w = c(1, 1, 0.8),
    warm_start = TRUE
  )
  
  fit_s <- fit_sparse(
    design_train, bases, ginfo, fit_p, M_grp,
    lambda_beta = tuning$lam_sp_b,
    lambda_alpha = tuning$lam_sp_a,
    lambda_group = tuning$lam_sp_g,
    ridge_beta = 2e-5,
    ridge_alpha = 1e-5,
    eta_adapt = tuning$eta_adapt,
    rough_w = c(1, 1, 0.8)
  )
  
  fit_r <- fit_refit(
    design_train, bases, ginfo, fit_s,
    lambda_beta = tuning$lam_rf_b,
    lambda_alpha = tuning$lam_rf_a,
    ridge_beta = 3e-5,
    ridge_alpha = 1e-5,
    active_tol = tuning$active_tol,
    active_method = tuning$active_method,
    min_gap = tuning$min_gap,
    slice_mix = 0.35,
    rough_w = c(1, 1, 0.8)
  )
  
  list(fit_p = fit_p, fit_s = fit_s, fit_r = fit_r, tuning = tuning)
}