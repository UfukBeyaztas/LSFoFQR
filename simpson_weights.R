simpson_weights <- function(x_grid) {
  x_grid <- validate_grid(x_grid, name = "x_grid")
  n <- length(x_grid)
  
  if (n %% 2 == 0) {
    stop("Simpson's rule requires an odd number of grid points. ",
         "Use integration = 'trapezoid' for general grids.")
  }
  
  dx <- diff(x_grid)
  if (max(abs(dx - dx[1])) > 1e-10) {
    stop("Simpson's rule requires an equally spaced grid. ",
         "Use integration = 'trapezoid' for non-equally spaced grids.")
  }
  
  h <- dx[1L]
  w <- c(1, rep(c(4, 2), (n - 3) / 2), 4, 1)
  w * h / 3
}