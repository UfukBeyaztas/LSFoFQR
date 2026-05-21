build_design_sparse <- function(Xmat, Ymat, grids, bases,
                                integration = c("trapezoid", "simpson", "auto")) {
  
  integration <- match.arg(integration)
  
  t_grid   <- grids$t
  s_grid   <- grids$s
  tau_grid <- grids$tau
  
  n  <- nrow(Xmat)
  ns <- length(s_grid)
  nu <- length(tau_grid)
  
  Kt  <- ncol(bases$Bt)
  Ks  <- ncol(bases$Bs)
  Ku  <- ncol(bases$Bu)
  Kas <- ncol(bases$Bsa)
  Kau <- ncol(bases$Bua)
  
  p_beta  <- Kt * Ks * Ku
  p_alpha <- Kas * Kau
  N       <- n * ns * nu
  
  R_scores <- compute_score_matrix(
    Xmat = Xmat,
    Bt = bases$Bt,
    t_grid = t_grid,
    integration = integration
  )
  
  active_s  <- lapply(seq_len(ns), function(l) which(bases$Bs[l, ]  != 0))
  active_u  <- lapply(seq_len(nu), function(r) which(bases$Bu[r, ]  != 0))
  active_sa <- lapply(seq_len(ns), function(l) which(bases$Bsa[l, ] != 0))
  active_ua <- lapply(seq_len(nu), function(r) which(bases$Bua[r, ] != 0))
  
  max_active_s  <- max(lengths(active_s))
  max_active_u  <- max(lengths(active_u))
  max_active_sa <- max(lengths(active_sa))
  max_active_ua <- max(lengths(active_ua))
  
  max_nnz_beta  <- as.numeric(N) * Kt * max_active_s * max_active_u
  max_nnz_alpha <- as.numeric(N) * max_active_sa * max_active_ua

  if (max_nnz_beta > .Machine$integer.max || max_nnz_alpha > .Machine$integer.max) {
    stop("Design too large for preallocated integer vectors. Reduce grid size or n.")
  }
  
  Xb_i <- integer(max_nnz_beta)
  Xb_j <- integer(max_nnz_beta)
  Xb_x <- numeric(max_nnz_beta)
  
  Xa_i <- integer(max_nnz_alpha)
  Xa_j <- integer(max_nnz_alpha)
  Xa_x <- numeric(max_nnz_alpha)
  
  y_vec   <- numeric(N)
  tau_vec <- numeric(N)
  ptr_b   <- 1
  ptr_a   <- 1
  row_id  <- 1
  
  for (r in seq_len(nu)) {
    bu_r  <- bases$Bu[r, ]
    bua_r <- bases$Bua[r, ]
    au_r  <- active_u[[r]]
    aua_r <- active_ua[[r]]
    
    for (l in seq_len(ns)) {
      bs_l  <- bases$Bs[l, ]
      bsa_l <- bases$Bsa[l, ]
      as_l  <- active_s[[l]]
      asa_l <- active_sa[[l]]
      
      col_alpha <- as.integer(outer(asa_l, aua_r, function(s2, u2) s2 + (u2 - 1) * Kas))
      val_alpha <- as.numeric(outer(bsa_l[asa_l], bua_r[aua_r]))
      n_a <- length(col_alpha)
      
      for (ii in seq_len(n)) {
        ri <- row_id + ii - 1
        end_a <- ptr_a + n_a - 1
        Xa_i[ptr_a:end_a] <- ri
        Xa_j[ptr_a:end_a] <- col_alpha
        Xa_x[ptr_a:end_a] <- val_alpha
        ptr_a <- end_a + 1
      }
      
      base_cols <- as.integer(outer(as_l, au_r,
                                    function(b, c) (b - 1) * Kt + (c - 1) * Kt * Ks))
      wsu <- as.numeric(outer(bs_l[as_l], bu_r[au_r]))
      n_su <- length(base_cols)
      n_b  <- Kt * n_su
      cols_pattern <- as.integer(rep(base_cols, each = Kt) + rep(seq_len(Kt), n_su))
      
      for (ii in seq_len(n)) {
        ri <- row_id + ii - 1
        end_b <- ptr_b + n_b - 1
        Xb_i[ptr_b:end_b] <- ri
        Xb_j[ptr_b:end_b] <- cols_pattern
        Xb_x[ptr_b:end_b] <- as.numeric(rep(wsu, each = Kt) * rep(R_scores[ii, ], n_su))
        ptr_b <- end_b + 1
      }
      
      rows_blk <- row_id:(row_id + n - 1)
      y_vec[rows_blk]   <- Ymat[, l]
      tau_vec[rows_blk] <- tau_grid[r]
      row_id <- row_id + n
    }
  }
  
  X_beta <- sparseMatrix(
    i = Xb_i[seq_len(ptr_b - 1)],
    j = Xb_j[seq_len(ptr_b - 1)],
    x = Xb_x[seq_len(ptr_b - 1)],
    dims = c(N, p_beta), repr = "C"
  )
  
  X_alpha <- sparseMatrix(
    i = Xa_i[seq_len(ptr_a - 1)],
    j = Xa_j[seq_len(ptr_a - 1)],
    x = Xa_x[seq_len(ptr_a - 1)],
    dims = c(N, p_alpha), repr = "C"
  )
  
  list(
    X_beta  = X_beta,
    X_alpha = X_alpha,
    y       = y_vec,
    tau_vec = tau_vec,
    R_scores = R_scores,
    N       = N,
    p_beta  = p_beta,
    p_alpha = p_alpha
  )
}
