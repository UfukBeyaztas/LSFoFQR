build_bspline_basis <- function(x, df, degree = 3, intercept = TRUE) {
  bs(x, df = df, degree = degree, intercept = intercept,
     Boundary.knots = range(x))
}
