generate_sfqfofr_data <- function(n_train = 100, n_test = 1000,
                                  t_grid, s_grid, tau_grid, bases,
                                  error_type = c("gaussian", "t3", "cauchy"),
                                  support_design = c("nonsparse", "nonsparse_noncrossing", "sparse"),
                                  integration = c("trapezoid", "simpson", "auto")) {
  
  integration <- match.arg(integration)
  error_type     <- match.arg(error_type)
  support_design <- match.arg(support_design)
  
  n_total <- n_train + n_test
  nt <- length(t_grid)
  ns <- length(s_grid)
  nu <- length(tau_grid)
  
  Kt  <- ncol(bases$Bt)
  Ks  <- ncol(bases$Bs)
  Ku  <- ncol(bases$Bu)
  Kas <- ncol(bases$Bsa)
  Kau <- ncol(bases$Bua)
  
  tau_bounded <- pmin(pmax(tau_grid, 1e-4), 1 - 1e-4)
  q_tau <- switch(
    error_type,
    gaussian = qnorm(tau_bounded),
    t3       = qt(tau_bounded, df = 3),
    cauchy   = qt(tau_bounded, df = 1)
  )
  
  c_const <- as.numeric(qr.solve(bases$Bu, rep(1, nu)))
  c_q     <- as.numeric(qr.solve(bases$Bu, q_tau))
  
  Theta_mean  <- matrix(0, Kt, Ks)
  Theta_scale <- matrix(0, Kt, Ks)
  
  patch <- function(mat, ir, jr, vals) {
    ir <- ir[ir >= 1 & ir <= nrow(mat)]
    jr <- jr[jr >= 1 & jr <= ncol(mat)]
    vals <- vals[seq_along(ir), seq_along(jr), drop = FALSE]
    mat[ir, jr] <- mat[ir, jr] + vals
    mat
  }
  
  if (support_design == "nonsparse") {
    if (Kt < 7 || Ks < 7) {
      stop("The nonsparse DGP needs Kt >= 7 and Ks >= 7. Recommended: kt=ks=11.")
    }
    
    ii <- seq_len(Kt)
    jj <- seq_len(Ks)
    
    damp <- outer(ii, jj, function(i, j) 1 / (1 + 0.20 * (i - 1) + 0.20 * (j - 1)))
    pat1 <- outer(ii, jj, function(i, j) sin(pi * i / (Kt + 1)) * cos(pi * j / (Ks + 1)))
    pat2 <- outer(ii, jj, function(i, j) cos(2 * pi * i / (Kt + 1)) * sin(pi * j / (Ks + 1)))
    pat3 <- outer(ii, jj, function(i, j) sin(2 * pi * i / (Kt + 1)) * sin(2 * pi * j / (Ks + 1)))
    pat4 <- outer(ii, jj, function(i, j) cos(pi * i / (Kt + 1)) * cos(2 * pi * j / (Ks + 1)))
    
    Theta_mean <- damp * (0.95 * pat1 + 0.65 * pat2 + 0.45 * pat3 + 0.30 * pat4)
    Theta_mean <- 1.70 * Theta_mean / max(abs(Theta_mean))
    
    Theta_scale <- damp * (0.55 * pat4 - 0.35 * pat2 + 0.25 * pat1)
    Theta_scale <- 0.10 * Theta_scale / max(abs(Theta_scale))
    
  } else if (support_design == "nonsparse_noncrossing") {
    if (Kt < 7 || Ks < 7) {
      stop("The nonsparse_noncrossing DGP needs Kt >= 7 and Ks >= 7. Recommended: kt=ks=11.")
    }
    
    ii <- seq_len(Kt)
    jj <- seq_len(Ks)
    
    damp <- outer(ii, jj, function(i, j) 1 / (1 + 0.12 * (i - 1) + 0.12 * (j - 1)))
    pos1 <- outer(ii, jj, function(i, j) 0.60 + 0.25 * sin(pi * i / (Kt + 1))^2 +
                    0.20 * cos(pi * j / (Ks + 1))^2)
    pos2 <- outer(ii, jj, function(i, j) 0.12 * sin(2 * pi * i / (Kt + 1)) *
                    sin(2 * pi * j / (Ks + 1)))
    Theta_mean <- damp * (pos1 + pos2)
    Theta_mean <- 1.35 * Theta_mean / mean(Theta_mean)
    
    raw_scale_pattern <- damp * outer(ii, jj, function(i, j) {
      0.50 * sin(pi * i / (Kt + 1)) * cos(2 * pi * j / (Ks + 1)) +
        0.35 * cos(2 * pi * i / (Kt + 1)) * sin(pi * j / (Ks + 1))
    })
    raw_scale_pattern <- raw_scale_pattern / max(abs(raw_scale_pattern))
    
    qmax <- max(abs(q_tau))
    max_allow <- 0.30 * min(Theta_mean) / max(qmax, 1e-8)
    Theta_scale <- min(0.055, max_allow) * raw_scale_pattern
    
    min_theta_over_tau <- min(vapply(q_tau, function(qq) min(Theta_mean + qq * Theta_scale), numeric(1L)))
    if (min_theta_over_tau <= 0) {
      stop("Internal DGP error: nonsparse_noncrossing coefficient matrix crosses zero.")
    }
    
  } else if (support_design == "sparse") {
    if (Kt < 11 || Ks < 11) {
      stop("The sparse DGP needs Kt >= 11 and Ks >= 11. Use kt=ks=11 or larger.")
    }
    
    idx_left_t  <- 2:4
    idx_left_s  <- 2:4
    idx_right_t <- (Kt - 3):(Kt - 1)
    idx_right_s <- (Ks - 3):(Ks - 1)
    
    bump_pos <- matrix(c(1.10, 1.45, 1.10,
                         1.45, 2.10, 1.45,
                         1.10, 1.45, 1.10), 3, 3, byrow = TRUE)
    bump_neg <- matrix(c(0.95, 1.30, 0.95,
                         1.30, 1.95, 1.30,
                         0.95, 1.30, 0.95), 3, 3, byrow = TRUE)
    
    Theta_mean <- patch(Theta_mean, idx_left_t,  idx_left_s,  bump_pos)
    Theta_mean <- patch(Theta_mean, idx_right_t, idx_right_s, -bump_neg)
    
    Theta_scale <- 0.08 * Theta_mean
  }
  
  a_mean <- rep(0, Kas)
  if (Kas >= 5) {
    a_mean[2] <- 0.08
    a_mean[4] <- -0.05
  } else if (Kas >= 2) {
    a_mean[2] <- 0.08
  } else {
    a_mean[1] <- 0.08
  }
  a_scale <- rep(0.25, Kas)
  
  phi_list <- list(
    sqrt(2) * sin(1 * pi * t_grid),
    sqrt(2) * cos(1 * pi * t_grid),
    sqrt(2) * sin(2 * pi * t_grid),
    sqrt(2) * cos(2 * pi * t_grid),
    sqrt(2) * sin(3 * pi * t_grid),
    sqrt(2) * cos(3 * pi * t_grid)
  )
  phi <- do.call(cbind, phi_list)
  score_sd <- c(2.20, 1.75, 1.35, 1.00, 0.75, 0.55)
  scores <- matrix(rnorm(n_total * length(score_sd)), n_total, length(score_sd))
  scores <- sweep(scores, 2, score_sd, `*`)
  
  Xmat <- scores %*% t(phi) + matrix(rnorm(n_total * nt, sd = 0.015), n_total, nt)
  Xmat <- scale(Xmat, center = TRUE, scale = FALSE)
  Xint <- compute_score_matrix(
    Xmat = Xmat,
    Bt = bases$Bt,
    t_grid = t_grid,
    integration = integration
  )
  
  alpha_mean_curve  <- as.numeric(bases$Bsa %*% a_mean)
  alpha_scale_curve <- as.numeric(bases$Bsa %*% a_scale)
  
  raw_mean  <- Xint %*% Theta_mean  %*% t(bases$Bs)
  raw_scale <- Xint %*% Theta_scale %*% t(bases$Bs)
  
  max_raw_scale   <- max(abs(raw_scale))
  min_alpha_scale <- min(alpha_scale_curve)
  x_scale_fac <- if (max_raw_scale > 0) {
    min(1, 0.75 * min_alpha_scale / max_raw_scale)
  } else {
    1
  }
  
  Xmat      <- x_scale_fac * Xmat
  Xint      <- x_scale_fac * Xint
  raw_mean  <- x_scale_fac * raw_mean
  raw_scale <- x_scale_fac * raw_scale
  
  mean_curve  <- raw_mean  + matrix(alpha_mean_curve,  n_total, ns, byrow = TRUE)
  scale_curve <- raw_scale + matrix(alpha_scale_curve, n_total, ns, byrow = TRUE)
  
  if (min(scale_curve) <= 0) stop("Scale curve not strictly positive; adjust DGP.")
  
  theta_true_arr <- array(0, c(Kt, Ks, Ku))
  for (k in seq_len(Ku)) {
    theta_true_arr[, , k] <- c_const[k] * Theta_mean + c_q[k] * Theta_scale
  }
  
  A_target   <- outer(a_mean, rep(1, nu)) + outer(a_scale, q_tau)
  Alpha_coef <- t(apply(A_target, 1, function(v) qr.solve(bases$Bua, v)))
  
  Sigma_s <- exp(-abs(outer(s_grid, s_grid, `-`)) / 0.24)
  diag(Sigma_s) <- 1
  
  if (error_type == "gaussian") {
    Eps <- mvrnorm(n_total, mu = rep(0, ns), Sigma = Sigma_s)
  } else {
    df_t <- if (error_type == "t3") 3 else 1
    Z <- mvrnorm(n_total, mu = rep(0, ns), Sigma = Sigma_s)
    scale_mix <- sqrt(df_t / rchisq(n_total, df = df_t))
    Eps <- Z * scale_mix
  }
  
  Ymat <- mean_curve + scale_curve * Eps
  
  beta_true_arr  <- array(0, c(nt, ns, nu))
  alpha_true_arr <- matrix(0, ns, nu)
  Qtrue          <- array(0, c(n_total, ns, nu))
  
  for (r in seq_len(nu)) {
    Theta_r <- Theta_mean + q_tau[r] * Theta_scale
    alpha_r <- alpha_mean_curve + q_tau[r] * alpha_scale_curve
    beta_true_arr[, , r] <- as.matrix(bases$Bt %*% Theta_r %*% t(bases$Bs))
    alpha_true_arr[, r] <- alpha_r
    Qtrue[, , r] <- mean_curve + scale_curve * q_tau[r]
  }
  
  tr <- seq_len(n_train)
  te <- (n_train + 1):n_total
  
  list(
    train = list(X = Xmat[tr, , drop = FALSE],
                 Y = Ymat[tr, , drop = FALSE],
                 Qtrue = Qtrue[tr, , , drop = FALSE]),
    test  = list(X = Xmat[te, , drop = FALSE],
                 Y = Ymat[te, , drop = FALSE],
                 Qtrue = Qtrue[te, , , drop = FALSE]),
    truth = list(
      Theta_mean = Theta_mean,
      Theta_scale = Theta_scale,
      theta_true_arr = theta_true_arr,
      alpha_coef = Alpha_coef,
      beta_arr = beta_true_arr,
      alpha_arr = alpha_true_arr,
      mean_curve = mean_curve,
      scale_curve = scale_curve,
      q_tau = q_tau,
      support_design = support_design,
      error_type = error_type
    )
  )
}