roughness_penalty_2d <- function(k1, k2, l1 = 1, l2 = 1) {
  P1 <- second_diff_penalty(k1)
  P2 <- second_diff_penalty(k2)
  l1 * kronecker(diag(k2), P1) + l2 * kronecker(P2, diag(k1))
}