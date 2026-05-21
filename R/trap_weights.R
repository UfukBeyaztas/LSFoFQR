trap_weights <- function(x_grid) {
  x_grid <- validate_grid(x_grid, name = "x_grid")
  n <- length(x_grid)
  
  w <- numeric(n)
  dx <- diff(x_grid)
  
  w[1] <- dx[1] / 2
  w[n]  <- dx[n - 1] / 2
  
  if (n > 2L) {
    w[2:(n - 1)] <- (dx[1:(n - 2)] + dx[2:(n - 1)]) / 2
  }
  
  w
}
