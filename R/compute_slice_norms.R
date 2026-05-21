compute_slice_norms <- function(beta_vec, ginfo) {
  G  <- ginfo$G
  Ku <- length(ginfo$slice_membership[[1]])
  balance <- ginfo$balance
  out <- matrix(0.0, G, Ku)
  
  for (g in seq_len(G)) {
    sm <- ginfo$slice_membership[[g]]
    for (k in seq_len(Ku)) {
      idx <- sm[[k]]
      if (length(idx) > 0) {
        out[g, k] <- sqrt(sum((balance[idx] * beta_vec[idx])^2))
      }
    }
  }
  out
}
