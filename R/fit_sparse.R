fit_sparse <- function(design, bases, ginfo, pilot_fit, M_grp,
                       lambda_beta = 1.5e-4, lambda_alpha = 6e-5,
                       lambda_group = 4e-3,
                       ridge_beta = 2e-5, ridge_alpha = 1e-5,
                       eta_adapt = 1.25, eps_adapt = 1e-4,
                       weight_cap = 1e4,
                       rough_w = c(1, 1, 0.8), solver = "L-BFGS-B",
                       eps_loss = 1e-6,
                       eps_grp = 1e-4,
                       maxit = 2000, factr = 1e7) {
  
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
  
  init_b <- pilot_fit$beta
  bn_pilot <- vapply(ginfo$groups, function(idx) {
    sqrt(sum((ginfo$balance[idx] * init_b[idx])^2))
  }, numeric(1))
  
  w_raw <- 1 / (bn_pilot + eps_adapt)^eta_adapt
  w_raw <- pmin(w_raw, weight_cap)
  w_g   <- lambda_group * w_raw
  
  fn <- function(par) {
    b <- par[seq_len(p_beta)]
    a <- par[p_beta + seq_len(p_alpha)]
    L <- smoothed_pinball_grad(b, a, design, eps_loss)
    pen <- as.numeric(t(b) %*% Q_b %*% b) + as.numeric(t(a) %*% Q_a %*% a)
    grp_sq  <- as.numeric(M_grp %*% (b^2))
    grp_pen <- sum(w_g * sqrt(grp_sq + eps_grp))
    L$loss + pen + grp_pen
  }
  
  gr <- function(par) {
    b <- par[seq_len(p_beta)]
    a <- par[p_beta + seq_len(p_alpha)]
    L <- smoothed_pinball_grad(b, a, design, eps_loss)
    
    gb <- -as.numeric(crossprod(design$X_beta,  L$grad_e)) + 2 * as.numeric(Q_b %*% b)
    ga <- -as.numeric(crossprod(design$X_alpha, L$grad_e)) + 2 * as.numeric(Q_a %*% a)
    
    grp_sq <- as.numeric(M_grp %*% (b^2))
    grad_grp_factor <- w_g / sqrt(grp_sq + eps_grp)
    gb_grp <- b * as.numeric(crossprod(M_grp, grad_grp_factor))
    
    c(gb + gb_grp, ga)
  }
  
  init <- c(pilot_fit$beta, pilot_fit$alpha)
  
  res <- optim(init, fn, gr, method = solver,
               control = list(maxit = maxit, factr = factr))
  
  b_hat <- res$par[seq_len(p_beta)]
  
  block_norms <- vapply(ginfo$groups, function(idx) {
    sqrt(sum((ginfo$balance[idx] * b_hat[idx])^2))
  }, numeric(1L))
  
  list(
    beta = b_hat,
    alpha = res$par[p_beta + seq_len(p_alpha)],
    status = if (res$convergence == 0) "optimal" else "converged_with_warning",
    value = res$value,
    convergence = res$convergence,
    message = res$message,
    block_norms = block_norms,
    pilot_block_norms = bn_pilot,
    adaptive_weights = w_raw,
    n_eval = res$counts[["function"]]
  )
}
