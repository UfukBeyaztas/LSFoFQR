default_tuning <- function(n_train) {

  if (missing(n_train) || length(n_train) != 1L ||
      !is.finite(n_train) || n_train <= 0) {
    stop("n_train must be a positive numeric value.")
  }

  n_train <- as.numeric(n_train)

  smooth_scale <- (100 / n_train)^0.25
  smooth_scale <- min(1.25, max(0.55, smooth_scale))

  group_scale <- (n_train / 100)^0.35

  active_tol <- 0.20 + 0.055 * log(n_train / 100)
  active_tol <- min(0.40, max(0.10, active_tol))

  list(
    ## Pilot fit
    lam_init_b = 4.0e-4 * smooth_scale,
    lam_init_a = 1.0e-4 * smooth_scale,

    ## Sparse fit
    lam_sp_b   = 1.5e-4 * smooth_scale,
    lam_sp_a   = 6.0e-5 * smooth_scale,
    lam_sp_g   = 1.0e-2 * group_scale,

    ## Constrained refit
    lam_rf_b   = 7.0e-5 * smooth_scale,
    lam_rf_a   = 5.0e-5 * smooth_scale,

    ## Active-set screening
    active_tol    = active_tol,
    eta_adapt     = if (n_train >= 250) 1.5 else 1.25,
    active_method = "fixed_or_gap",
    min_gap       = 0.02
  )
}
