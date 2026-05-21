select_blocks_by_score <- function(scores,
                                   method = c("gap", "fixed", "fixed_or_gap"),
                                   active_tol = 0.20,
                                   min_gap = 0.02,
                                   min_active = 1) {
  method <- match.arg(method)
  G <- length(scores)
  
  if (G == 1) {
    return(list(active_blocks = 1, kappa = -Inf, method_used = "single"))
  }
  
  fixed_sel <- which(scores > active_tol)
  if (length(fixed_sel) < min_active) fixed_sel <- order(scores, decreasing = TRUE)[seq_len(min_active)]
  
  ord <- order(scores, decreasing = TRUE)
  scr <- scores[ord]
  gaps <- scr[-G] - scr[-1]
  k_gap <- which.max(gaps)
  gap_size <- gaps[k_gap]
  
  if (is.na(gap_size) || gap_size < min_gap) {
    gap_sel <- fixed_sel
    kappa_gap <- active_tol
    gap_used <- FALSE
  } else {
    kappa_gap <- 0.5 * (scr[k_gap] + scr[k_gap + 1])
    gap_sel <- which(scores > kappa_gap)
    if (length(gap_sel) < min_active) gap_sel <- ord[seq_len(min_active)]
    gap_used <- TRUE
  }
  
  if (method == "fixed") {
    active <- fixed_sel
    kappa  <- active_tol
    used   <- "fixed"
  } else if (method == "gap") {
    active <- gap_sel
    kappa  <- kappa_gap
    used   <- if (gap_used) "gap" else "fixed_fallback"
  } else {
    active <- union(fixed_sel, gap_sel)
    kappa  <- min(active_tol, kappa_gap)
    used   <- if (gap_used) "fixed_or_gap" else "fixed_fallback"
  }
  
  list(active_blocks = sort(active), kappa = kappa, method_used = used,
       gap_size = gap_size, gap_rank = k_gap)
}
