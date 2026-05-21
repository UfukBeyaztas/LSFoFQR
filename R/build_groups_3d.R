build_groups_3d <- function(d1, d2, d3, b1 = 3, b2 = 3, b3 = 3) {
  n1 <- d1 - b1 + 1
  n2 <- d2 - b2 + 1
  n3 <- d3 - b3 + 1
  groups <- vector("list", n1 * n2 * n3)
  g <- 1
  
  for (k in seq_len(n3)) {
    for (j in seq_len(n2)) {
      for (i in seq_len(n1)) {
        ii <- i:(i + b1 - 1)
        jj <- j:(j + b2 - 1)
        kk <- k:(k + b3 - 1)
        groups[[g]] <- as.integer(
          outer(ii,
                as.integer(outer(jj, kk,
                                 function(jv, kv) (jv - 1) * d1 + (kv - 1) * d1 * d2)),
                `+`)
        )
        g <- g + 1
      }
    }
  }
  groups
}
