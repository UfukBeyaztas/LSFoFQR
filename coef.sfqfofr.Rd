\name{coef.sfqfofr}
\alias{coef.sfqfofr}
\alias{coef}
\title{Extract Coefficient Surfaces from an \code{sfqfofr} Fit}

\description{
Extracts the estimated coefficient surface from a fitted simultaneous
function-on-function quantile regression model of class \code{"sfqfofr"}.
The returned object is the estimated regression kernel
\eqn{\widehat{\beta}(t,s,u)} evaluated on the predictor grid, response grid,
and quantile grid used in the fitted model.
}

\usage{
\method{coef}{sfqfofr}(object, tau_index = NULL, ...)
}

\arguments{
  \item{object}{
    A fitted model object of class \code{"sfqfofr"}, typically returned by
    \code{\link{sfqfofr}}. The object must contain the refitted coefficient
    vector, the B-spline basis objects, and the grids used for estimation.
  }

  \item{tau_index}{
    Optional integer index specifying the quantile level for which the
    coefficient surface should be returned. If \code{NULL}, the full
    three-dimensional array of coefficient surfaces over all quantile levels is
    returned. If supplied, \code{tau_index} should be an integer between
    \code{1} and \code{length(object$grids$tau)}.
  }

  \item{...}{
    Additional arguments passed to the S3 method. Currently unused.
  }
}

\details{
The fitted simultaneous function-on-function quantile regression model has the
form
\deqn{
Q_{Y(s)\mid X}(u \mid X)
=
\alpha(s,u)
+
\int_{\mathcal{T}} X(t)\beta(t,s,u)\,dt,
\qquad (s,u)\in\mathcal{S}\times\mathcal{U},
}
where \eqn{X(t)} is the functional predictor, \eqn{Y(s)} is the functional
response, \eqn{\alpha(s,u)} is the quantile-specific intercept surface, and
\eqn{\beta(t,s,u)} is the quantile-varying function-on-function regression
kernel. The quantity \eqn{\beta(t,s,u)} measures the contribution of the
predictor trajectory at location \eqn{t} to the \eqn{u}th conditional quantile
of the response at location \eqn{s}.

In the fitted model, the coefficient kernel is approximated by a tensor-product
B-spline expansion,
\deqn{
\beta(t,s,u)
\approx
\sum_{a=1}^{K_t}\sum_{b=1}^{K_s}\sum_{c=1}^{K_u}
\theta_{abc} b_a(t)c_b(s)d_c(u),
}
where \eqn{\{b_a(t)\}}, \eqn{\{c_b(s)\}}, and \eqn{\{d_c(u)\}} are B-spline
bases over the predictor, response, and quantile domains, respectively. The
coefficient tensor \eqn{\theta_{abc}} is estimated by the three-stage procedure
implemented in \code{\link{sfqfofr}}: a smooth pilot fit, an adaptive locally
sparse fit, and a constrained refit on the selected active coefficient set.

This function reconstructs the estimated coefficient surface
\eqn{\widehat{\beta}(t,s,u)} from the final refitted coefficient vector. If
\code{tau_index = NULL}, the function evaluates the fitted kernel at all
quantile levels in \code{object$grids$tau}. If \code{tau_index} is supplied,
only the corresponding quantile-specific coefficient surface is returned.

The returned coefficient surface is evaluated on the same predictor and response
grids used in the fitted object. It is therefore directly suitable for plotting,
visual inspection of locally inactive regions, and interpretation of the
lag-response-quantile dependence structure. In particular, regions where the
estimated surface is exactly zero correspond to regions selected as locally
inactive by the overlapping group-Lasso component of the estimator.
}

\value{
If \code{tau_index = NULL}, a three-dimensional numeric array with dimensions
\deqn{
\texttt{length(object\$grids\$t)}
\times
\texttt{length(object\$grids\$s)}
\times
\texttt{length(object\$grids\$tau)}.
}
The first dimension corresponds to the predictor-domain grid, the second to the
response-domain grid, and the third to the quantile grid.

If \code{tau_index} is supplied, the function returns a three-dimensional array
with the same first two dimensions and a singleton third dimension. Thus, the
returned object has dimension
\deqn{
\texttt{length(object\$grids\$t)}
\times
\texttt{length(object\$grids\$s)}
\times
1.
}
}
\author{
Ufuk Beyaztas
}

\note{
This function extracts only the estimated coefficient kernel
\eqn{\widehat{\beta}(t,s,u)}. To obtain predicted conditional quantile response
curves for new functional predictors, use \code{\link{predict.sfqfofr}}.
}
\examples{
\donttest{
## Grids and basis

t_grid   <- seq(0, 1, length.out = 101)
s_grid   <- seq(0, 1, length.out = 101)
tau_grid <- c(0.10, 0.35, 0.50, 0.65, 0.80, 0.90)

bases <- make_sfqfofr_basis(
  t_grid = t_grid,
  s_grid = s_grid,
  tau_grid = tau_grid,
  kt = 11,
  ks = 11,
  ku = length(tau_grid),
  kas = 5,
  kau = length(tau_grid),
  degree_t = 2,
  degree_s = 2,
  degree_u = 2
)

## Generate data

dat <- generate_sfqfofr_data(
  n_train = 100,
  n_test = 1000,
  t_grid = t_grid,
  s_grid = s_grid,
  tau_grid = tau_grid,
  bases = bases,
  error_type = "gaussian",
  support_design = "sparse",
  integration = "trapezoid"
)

## Sparse-design-friendly tuning values

tuning_sparse_example <- list(
  lam_init_b = 2.0e-04,
  lam_init_a = 5.0e-05,

  lam_sp_b   = 7.5e-05,
  lam_sp_a   = 3.0e-05,
  lam_sp_g   = 2.5e-02,

  lam_rf_b   = 3.5e-05,
  lam_rf_a   = 2.5e-05,

  active_tol    = 0.32,
  eta_adapt     = 1.5,
  active_method = "fixed_or_gap",
  min_gap       = 0.02
)

## Fit sfqfofr model

fit <- sfqfofr(
  Xmat = dat$train$X,
  Ymat = dat$train$Y,
  grids = list(t = t_grid, s = s_grid, tau = tau_grid),
  bases = bases,
  tuning = FALSE,
  base_tuning = tuning_sparse_example,
  integration = "trapezoid"
)

## Extract true and estimated beta(t,s,u) at u = 0.50

tau_id <- which.min(abs(tau_grid - 0.50))

beta_true_med <- dat$truth$beta_arr[, , tau_id]
beta_hat_med  <- coef(fit, tau_index = tau_id)[, , 1]

## Common z-axis range for fair visual comparison

zlim_common <- range(c(beta_true_med, beta_hat_med), na.rm = TRUE)

zpad <- 0.05 * diff(zlim_common)
if (!is.finite(zpad) || zpad == 0) zpad <- 0.1

zlim_common <- zlim_common + c(-zpad, zpad)

## Plot true and estimated coefficient surfaces side by side

op <- par(
  mfrow = c(1, 2),
  mar = c(2.5, 2.5, 3.0, 1.0),
  oma = c(0.2, 0.2, 0.2, 0.2)
)

persp(
  x = t_grid,
  y = s_grid,
  z = beta_true_med,
  theta = 35,
  phi = 25,
  expand = 0.75,
  col = "gray92",
  border = "gray70",
  ticktype = "detailed",
  xlab = "t",
  ylab = "s",
  zlab = expression(beta(t, s)),
  zlim = zlim_common,
  main = expression(True~~beta(t, s)~~at~~u == 0.50)
)

persp(
  x = t_grid,
  y = s_grid,
  z = beta_hat_med,
  theta = 35,
  phi = 25,
  expand = 0.75,
  col = "gray92",
  border = "gray70",
  ticktype = "detailed",
  xlab = "t",
  ylab = "s",
  zlab = expression(hat(beta)(t, s)),
  zlim = zlim_common,
  main = expression(Estimated~~beta(t, s)~~at~~u == 0.50)
)

par(op)
}
}
