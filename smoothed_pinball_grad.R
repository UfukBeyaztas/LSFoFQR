smoothed_pinball_grad <- function(b, a, design, eps_loss = 1e-6) {
  e <- design$y - as.numeric(design$X_beta %*% b) - as.numeric(design$X_alpha %*% a)
  s <- sqrt(e^2 + eps_loss)
  loss <- sum(0.5 * s + (design$tau_vec - 0.5) * e) / design$N
  grad_e <- (0.5 * e / s + (design$tau_vec - 0.5)) / design$N
  list(loss = loss, grad_e = grad_e, residual = e)
}