empirical_check_loss_design <- function(fit_r, design) {
  e <- design$y - as.numeric(design$X_beta %*% fit_r$beta) -
    as.numeric(design$X_alpha %*% fit_r$alpha)
  mean(0.5 * abs(e) + (design$tau_vec - 0.5) * e)
}
