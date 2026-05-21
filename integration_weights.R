integration_weights <- function(x_grid,
                                method = c("trapezoid", "simpson", "auto")) {
  
  method <- match.arg(method)
  x_grid <- validate_grid(x_grid, name = "x_grid")
  
  if (method == "trapezoid") {
    return(trap_weights(x_grid))
  }
  
  if (method == "simpson") {
    return(simpson_weights(x_grid))
  }
  
  n <- length(x_grid)
  dx <- diff(x_grid)
  equally_spaced <- max(abs(dx - dx[1L])) <= 1e-10
  
  if (n %% 2 == 1 && equally_spaced) {
    return(simpson_weights(x_grid))
  }
  
  trap_weights(x_grid)
}