build_group_matrix <- function(ginfo, p_beta) {
  G <- ginfo$G
  i_idx <- rep(seq_len(G), lengths(ginfo$groups))
  j_idx <- unlist(ginfo$groups, use.names = FALSE)
  x_val <- (ginfo$balance[j_idx])^2
  sparseMatrix(i = i_idx, j = j_idx, x = x_val,
               dims = c(G, p_beta), repr = "C")
}
