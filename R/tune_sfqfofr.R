tune_sfqfofr <- function(design_train, bases, ginfo, M_grp,
                         n_subjects,
                         tuning_grid = NULL,
                         base_tuning = default_tuning(n_subjects),
                         criterion = "qbic",
                         bic_c = 1.0,
                         df_type = "active_coeff") {
  criterion <- match.arg(criterion, choices = "qbic")
  if (is.null(tuning_grid)) tuning_grid <- make_tuning_grid(n_subjects, base_tuning)
  
  results <- vector("list", length(tuning_grid))
  fits <- vector("list", length(tuning_grid))
  
  for (jj in seq_along(tuning_grid)) {
    tun <- tuning_grid[[jj]]
    
    fit_j <- tryCatch(
      sfqfofr_tuning(design_train, bases, ginfo, M_grp, tun),
      error = function(e) {
        NULL
      }
    )
    
    if (is.null(fit_j)) {
      results[[jj]] <- cbind(tuning_list_to_row(tun, jj),
                             criterion_value = Inf, train_check = NA_real_,
                             qbic = NA_real_,
                             active_blocks = NA_integer_, active_coeff = NA_integer_,
                             status = "failed")
      fits[[jj]] <- NULL
      next
    }
    
    qinfo <- qbic_sfqfofr(fit_j$fit_r, design_train, n_subjects,
                          bic_c = bic_c, df_type = df_type)
    crit_val <- qinfo$qbic
    
    results[[jj]] <- cbind(
      tuning_list_to_row(tun, jj),
      criterion_value = crit_val,
      train_check = qinfo$loss,
      qbic = qinfo$qbic,
      df_total = qinfo$df_total,
      active_blocks = length(fit_j$fit_r$active_blocks),
      active_coeff = length(fit_j$fit_r$active_idx),
      status = fit_j$fit_r$status
    )
    fits[[jj]] <- fit_j
  }
  
  tab <- do.call(rbind, results)
  ok <- which(is.finite(tab$criterion_value))
  if (!length(ok)) stop("All tuning candidates failed.")
  best_id <- ok[which.min(tab$criterion_value[ok])]
  
  list(
    best_id = best_id,
    best_tuning = tuning_grid[[best_id]],
    best_fit = fits[[best_id]],
    tuning_table = tab
  )
}
