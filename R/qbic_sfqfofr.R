qbic_sfqfofr <- function(fit_r, design, n_subjects,
                         bic_c = 1.0,
                         df_type = c("active_coeff", "active_blocks")) {
  df_type <- match.arg(df_type)
  loss <- empirical_check_loss_design(fit_r, design)
  df_beta <- if (df_type == "active_coeff") {
    length(fit_r$active_idx)
  } else {
    length(fit_r$active_blocks)
  }
  df_total <- df_beta + design$p_alpha
  qbic <- log(loss + 1e-12) + bic_c * df_total * log(n_subjects) / n_subjects
  list(qbic = qbic, loss = loss, df_total = df_total, df_beta = df_beta)
}
