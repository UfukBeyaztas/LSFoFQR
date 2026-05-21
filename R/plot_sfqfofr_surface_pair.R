plot_sfqfofr_surface_pair <- function(beta_true, beta_hat, t_grid, s_grid,
                                      tau_value = 0.50) {
  zlim_common <- range(c(beta_true, beta_hat), na.rm = TRUE)
  zpad <- 0.05 * diff(zlim_common)
  if (!is.finite(zpad) || zpad == 0) zpad <- 0.1
  zlim_common <- zlim_common + c(-zpad, zpad)

  op <- par(mfrow = c(1, 2), mar = c(2.5, 2.5, 2.5, 1.0))

  persp(
    x = t_grid, y = s_grid, z = beta_true,
    theta = 35, phi = 25, expand = 0.75,
    col = "gray92", border = "gray70",
    ticktype = "detailed",
    xlab = "t", ylab = "s", zlab = expression(beta(t, s)),
    zlim = zlim_common,
    main = bquote(True ~ beta(t, s) ~ "at" ~ u == .(tau_value))
  )

  persp(
    x = t_grid, y = s_grid, z = beta_hat,
    theta = 35, phi = 25, expand = 0.75,
    col = "gray92", border = "gray70",
    ticktype = "detailed",
    xlab = "t", ylab = "s", zlab = expression(beta(t, s)),
    zlim = zlim_common,
    main = bquote(Estimated ~ beta(t, s) ~ "at" ~ u == .(tau_value))
  )

  par(op)
}
