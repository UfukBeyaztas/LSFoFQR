fit_refit <- function(design, bases, ginfo, sparse_fit,
                      lambda_beta = 7e-5, lambda_alpha = 5e-5,
                      ridge_beta = 3e-5, ridge_alpha = 1e-5,
                      active_tol = 0.20, slice_mix = 0.35,
                      active_method = c("gap", "fixed", "fixed_or_gap"),
                      min_gap = 0.02,
                      rough_w = c(1, 1, 0.8), solver = "L-BFGS-B",
                      eps_loss = 1e-6,
                      maxit = 2000, factr = 1e7) {
  
  active_method <- match.arg(active_method)
  
  Kt  <- ncol(bases$Bt)
  Ks  <- ncol(bases$Bs)
  Ku  <- ncol(bases$Bu)
  Kas <- ncol(bases$Bsa)
  Kau <- ncol(bases$Bua)
  
  p_beta  <- design$p_beta
  p_alpha <- design$p_alpha
  
  groups <- ginfo$groups
  bn <- sparse_fit$block_norms
  
  G_g <- bn / max(bn, 1e-8)
  
  sn <- compute_slice_norms(sparse_fit$beta, ginfo)
  sm <- apply(sn, 2, max)
  sm[sm < 1e-8] <- 1.0
  A_g <- rowMeans(sweep(sn, 2, sm, "/"))
  
  H_g <- (1 - slice_mix) * G_g + slice_mix * A_g
  
  sel <- select_blocks_by_score(H_g, method = active_method,
                                active_tol = active_tol, min_gap = min_gap)
  active_blocks <- sel$active_blocks
  
  active_idx <- sort(unique(unlist(groups[active_blocks], use.names = FALSE)))
  inactive_idx <- setdiff(seq_len(p_beta), active_idx)
  p_act <- length(active_idx)
  
  X_beta_act <- design$X_beta[, active_idx, drop = FALSE]
  
  P_beta_full <- roughness_penalty_3d(Kt, Ks, Ku, rough_w[1], rough_w[2], rough_w[3])
  P_beta_act  <- P_beta_full[active_idx, active_idx, drop = FALSE]
  P_alpha     <- roughness_penalty_2d(Kas, Kau, 1, 0.5)
  
  Q_b_act <- lambda_beta  * P_beta_act + ridge_beta  * diag(p_act)
  Q_a     <- lambda_alpha * P_alpha    + ridge_alpha * diag(p_alpha)
  
  fn <- function(par) {
    b <- par[seq_len(p_act)]
    a <- par[p_act + seq_len(p_alpha)]
    e <- design$y - as.numeric(X_beta_act %*% b) - as.numeric(design$X_alpha %*% a)
    loss <- sum(0.5 * sqrt(e^2 + eps_loss) + (design$tau_vec - 0.5) * e) / design$N
    pen <- as.numeric(t(b) %*% Q_b_act %*% b) + as.numeric(t(a) %*% Q_a %*% a)
    loss + pen
  }
  
  gr <- function(par) {
    b <- par[seq_len(p_act)]
    a <- par[p_act + seq_len(p_alpha)]
    e <- design$y - as.numeric(X_beta_act %*% b) - as.numeric(design$X_alpha %*% a)
    grad_e <- (0.5 * e / sqrt(e^2 + eps_loss) + (design$tau_vec - 0.5)) / design$N
    gb <- -as.numeric(crossprod(X_beta_act, grad_e)) + 2 * as.numeric(Q_b_act %*% b)
    ga <- -as.numeric(crossprod(design$X_alpha, grad_e)) + 2 * as.numeric(Q_a %*% a)
    c(gb, ga)
  }
  
  init <- c(sparse_fit$beta[active_idx], sparse_fit$alpha)
  
  res <- optim(init, fn, gr, method = solver,
               control = list(maxit = maxit, factr = factr))
  
  beta_full <- numeric(p_beta)
  beta_full[active_idx] <- res$par[seq_len(p_act)]
  
  list(
    beta = beta_full,
    alpha = res$par[p_act + seq_len(p_alpha)],
    status = if (res$convergence == 0) "optimal" else "converged_with_warning",
    value = res$value,
    convergence = res$convergence,
    message = res$message,
    active_blocks = active_blocks,
    active_idx = active_idx,
    inactive_idx = inactive_idx,
    global_scores = G_g,
    slice_scores = A_g,
    combined_scores = H_g,
    kappa_used = sel$kappa,
    active_method_used = sel$method_used,
    gap_size = sel$gap_size,
    gap_rank = sel$gap_rank,
    n_eval = res$counts[["function"]]
  )
}