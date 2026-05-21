build_group_info <- function(groups, p, Kt, Ks, Ku) {
  mult    <- coefficient_multiplicity(groups, p)
  balance <- 1 / sqrt(mult)
  p_slice <- Kt * Ks
  
  slice_membership <- lapply(groups, function(idx) {
    kk <- ((idx - 1) %/% p_slice) + 1
    lapply(seq_len(Ku), function(k) idx[kk == k])
  })
  
  list(
    groups = groups,
    mult = mult,
    balance = balance,
    slice_membership = slice_membership,
    G = length(groups)
  )
}
