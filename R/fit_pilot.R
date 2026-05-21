fit_pilot <- function(design, bases, ginfo,
                      lambda_beta = 4e-4, lambda_alpha = 1e-4,
                      ridge_beta = 2e-5, ridge_alpha = 1e-5,
                      rough_w = c(1, 1, 0.8), solver = "L-BFGS-B",
                      eps_loss = 1e-6,
                      warm_start = TRUE,
                      maxit = 1000, factr = 1e7) {
  
  Kt  <- ncol(bases$Bt)
  Ks  <- ncol(bases$Bs)
  Ku  <- ncol(bases$Bu)
  Kas <- ncol(bases$Bsa)
  Kau <- ncol(bases$Bua)
  
  p_beta  <- design$p_beta
  p_alpha <- design$p_alpha
  
  P_beta  <- roughness_penalty_3d(Kt, Ks, Ku, rough_w[1], rough_w[2], rough_w[3])
  P_alpha <- roughness_penalty_2d(Kas, Kau, 1, 0.5)
  
  Q_b <- lambda_beta  * P_beta  + ridge_beta  * diag(p_beta)
  Q_a <- lambda_alpha * P_alpha + ridge_alpha * diag(p_alpha)
  
  fn <- function(par) {
    b <- par[seq_len(p_beta)]
    a <- par[p_beta + seq_len(p_alpha)]
    L <- smoothed_pinball_grad(b, a, design, eps_loss)
    L$loss + as.numeric(t(b) %*% Q_b %*% b) + as.numeric(t(a) %*% Q_a %*% a)
  }
  
  gr <- function(par) {
    b <- par[seq_len(p_beta)]
    a <- par[p_beta + seq_len(p_alpha)]
    L <- smoothed_pinball_grad(b, a, design, eps_loss)
    gb <- -as.numeric(crossprod(design$X_beta,  L$grad_e)) + 2 * as.numeric(Q_b %*% b)
    ga <- -as.numeric(crossprod(design$X_alpha, L$grad_e)) + 2 * as.numeric(Q_a %*% a)
    c(gb, ga)
  }
  
  init <- if (warm_start) {
    make_ridge_warm_start(design, Q_b, Q_a)
  } else {
    rep(0, p_beta + p_alpha)
  }
  
  res <- optim(init, fn, gr, method = solver,
               control = list(maxit = maxit, factr = factr))
  
  list(
    beta = res$par[seq_len(p_beta)],
    alpha = res$par[p_beta + seq_len(p_alpha)],
    status = if (res$convergence == 0) "optimal" else "converged_with_warning",
    value = res$value,
    convergence = res$convergence,
    message = res$message,
    n_eval = res$counts[["function"]]
  )
}
