sfqfofr <- function(Xmat, Ymat, grids, bases,
                    tuning = FALSE,
                    tuning_grid = NULL,
                    tuning_criterion = "qbic",
                    base_tuning = NULL,
                    bic_c = 1.0,
                    df_type = "active_coeff",
                    integration = c("trapezoid", "simpson", "auto")) {
  
  integration <- match.arg(integration)
  tuning_criterion <- match.arg(tuning_criterion, choices = "qbic")
  n_train <- nrow(Xmat)
  if (is.null(base_tuning)) base_tuning <- default_tuning(n_train)
  
  if (ncol(Xmat) != length(grids$t)) {
    stop("ncol(Xmat) must equal length(grids$t).")
  }
  
  if (ncol(Ymat) != length(grids$s)) {
    stop("ncol(Ymat) must equal length(grids$s).")
  }
  
  if (length(grids$tau) < 2L) {
    stop("At least two quantile levels are recommended for simultaneous quantile estimation.")
  }
  
  design_train <- build_design_sparse(
    Xmat = Xmat,
    Ymat = Ymat,
    grids = grids,
    bases = bases,
    integration = integration
  )
  
  Kt <- ncol(bases$Bt); Ks <- ncol(bases$Bs); Ku <- ncol(bases$Bu)
  groups3d <- build_groups_3d(Kt, Ks, Ku, 3, 3, 3)
  ginfo <- build_group_info(groups3d, p = Kt * Ks * Ku, Kt = Kt, Ks = Ks, Ku = Ku)
  M_grp <- build_group_matrix(ginfo, p_beta = Kt * Ks * Ku)
  
  if (tuning) {
    tuned <- tune_sfqfofr(
      design_train = design_train,
      bases = bases,
      ginfo = ginfo,
      M_grp = M_grp,
      n_subjects = n_train,
      tuning_grid = tuning_grid,
      base_tuning = base_tuning,
      criterion = tuning_criterion,
      bic_c = bic_c,
      df_type = df_type
    )
    fit_object <- tuned$best_fit
    tuning_used <- tuned$best_tuning
    tuning_table <- tuned$tuning_table
  } else {
    fit_object <- sfqfofr_tuning(
      design_train = design_train,
      bases = bases,
      ginfo = ginfo,
      M_grp = M_grp,
      tuning = base_tuning
    )
    tuning_used <- base_tuning
    tuning_table <- NULL
  }
  
  out <- list(
    fit_p = fit_object$fit_p,
    fit_s = fit_object$fit_s,
    fit_r = fit_object$fit_r,
    tuning = tuning_used,
    tuning_table = tuning_table,
    design_train = design_train,
    ginfo = ginfo,
    M_grp = M_grp,
    bases = bases,
    grids = grids,
    integration = integration
  )
  
  class(out) <- "sfqfofr"
  
  out
}