validate_grid <- function(x_grid, name = "grid") {
  x_grid <- as.numeric(x_grid)
  
  if (length(x_grid) < 2) {
    stop(name, " must contain at least two grid points.")
  }
  
  if (any(!is.finite(x_grid))) {
    stop(name, " contains non-finite values.")
  }
  
  if (any(diff(x_grid) <= 0)) {
    stop(name, " must be strictly increasing.")
  }
  
  x_grid
}