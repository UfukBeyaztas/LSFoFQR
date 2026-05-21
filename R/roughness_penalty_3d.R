roughness_penalty_3d <- function(k1, k2, k3, l1 = 1, l2 = 1, l3 = 1) {
  P1 <- second_diff_penalty(k1)
  P2 <- second_diff_penalty(k2)
  P3 <- second_diff_penalty(k3)
  l1 * kronecker(diag(k3), kronecker(diag(k2), P1)) +
    l2 * kronecker(diag(k3), kronecker(P2, diag(k1))) +
    l3 * kronecker(P3, kronecker(diag(k2), diag(k1)))
}
