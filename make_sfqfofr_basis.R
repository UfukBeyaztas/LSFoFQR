make_sfqfofr_basis <- function(t_grid, s_grid, tau_grid,
                               kt = 11, ks = 11, ku = NULL,
                               kas = 5, kau = NULL,
                               degree_t = 2L, degree_s = 2, degree_u = 2) {
  if (is.null(ku))  ku  <- length(tau_grid)
  if (is.null(kau)) kau <- length(tau_grid)
  list(
    Bt  = build_bspline_basis(t_grid,   df = kt,  degree = degree_t, intercept = TRUE),
    Bs  = build_bspline_basis(s_grid,   df = ks,  degree = degree_s, intercept = TRUE),
    Bu  = build_bspline_basis(tau_grid, df = ku,  degree = degree_u, intercept = TRUE),
    Bsa = build_bspline_basis(s_grid,   df = kas, degree = degree_s, intercept = TRUE),
    Bua = build_bspline_basis(tau_grid, df = kau, degree = degree_u, intercept = TRUE)
  )
}
