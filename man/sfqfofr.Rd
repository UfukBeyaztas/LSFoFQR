\name{sfqfofr}
\alias{sfqfofr}
\title{Fit Locally Sparse Simultaneous Function-on-Function Quantile Regression}

\description{
Fits the proposed locally sparse simultaneous function-on-function quantile
regression model. The method estimates a quantile-varying regression kernel
\eqn{\beta(t,s,u)} jointly over the predictor domain, response domain, and
quantile domain using tensor-product B-splines, directional roughness penalties,
and an adaptive overlapping group-Lasso penalty.
}

\usage{
sfqfofr(
  Xmat,
  Ymat,
  grids,
  bases,
  tuning = FALSE,
  tuning_grid = NULL,
  tuning_criterion = "qbic",
  base_tuning = NULL,
  bic_c = 1.0,
  df_type = "active_coeff",
  integration = c("trapezoid", "simpson", "auto")
)
}

\arguments{
  \item{Xmat}{
    Numeric matrix of functional predictors. Each row corresponds to one
    functional predictor curve and each column corresponds to a grid point in
    \code{grids$t}. Thus, \code{ncol(Xmat)} must equal
    \code{length(grids$t)}.
  }

  \item{Ymat}{
    Numeric matrix of functional responses. Each row corresponds to one
    functional response curve and each column corresponds to a grid point in
    \code{grids$s}. Thus, \code{ncol(Ymat)} must equal
    \code{length(grids$s)}.
  }

  \item{grids}{
    A list containing the grids used for fitting. It must contain:
    \itemize{
      \item \code{t}: numeric vector giving the predictor-domain grid;
      \item \code{s}: numeric vector giving the response-domain grid;
      \item \code{tau}: numeric vector of quantile levels.
    }
    The quantile grid should contain at least two levels for simultaneous
    quantile estimation.
  }

  \item{bases}{
    A list of B-spline basis matrices, typically returned by
    \code{\link{make_sfqfofr_basis}}. The list must contain \code{Bt},
    \code{Bs}, \code{Bu}, \code{Bsa}, and \code{Bua}.
  }

  \item{tuning}{
    Logical value indicating whether tuning-parameter selection should be
    performed. If \code{FALSE}, the estimator is fitted using
    \code{base_tuning}. If \code{TRUE}, the candidate vectors in
    \code{tuning_grid} are evaluated and the best one is selected according to
    \code{tuning_criterion}.
  }

  \item{tuning_grid}{
    Optional list of candidate tuning-parameter vectors. Each element should be
    a named list containing the same components as \code{\link{default_tuning}}.
    If \code{NULL} and \code{tuning = TRUE}, a grid is constructed internally
    around \code{base_tuning}.
  }

  \item{tuning_criterion}{
    Character string specifying the tuning criterion. Currently only
    \code{"qbic"} is supported.
  }

  \item{base_tuning}{
    Named list of tuning parameters used as the baseline fitting vector. If
    \code{NULL}, \code{default_tuning(nrow(Xmat))} is used.
  }

  \item{bic_c}{
    Positive numeric multiplier used in the QBIC-type tuning criterion. Larger
    values impose a stronger penalty on model complexity.
  }

  \item{df_type}{
    Character string specifying how the effective degrees of freedom are counted
    in the QBIC-type criterion. The default is \code{"active_coeff"}, which
    counts the intercept coefficients and the selected active slope coefficients
    after support screening.
  }

  \item{integration}{
    Numerical integration rule used to approximate functional inner products.
    Available options are \code{"trapezoid"}, \code{"simpson"}, and
    \code{"auto"}. The default is \code{"trapezoid"}, which is recommended for
    general grids because it works for arbitrary grid lengths.
  }
}

\details{
The simultaneous function-on-function quantile regression model fitted by this
function is
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
\eqn{\beta(t,s,u)} is the quantile-varying function-on-function coefficient
kernel. The coefficient kernel \eqn{\beta(t,s,u)} describes how the predictor
value at location \eqn{t} affects the \eqn{u}th conditional quantile of the
response at location \eqn{s}.

The coefficient kernel is approximated by a tensor-product B-spline expansion,
\deqn{
\beta(t,s,u)
\approx
\sum_{a=1}^{K_t}
\sum_{b=1}^{K_s}
\sum_{c=1}^{K_u}
\theta_{abc} b_a(t)c_b(s)d_c(u),
}
where \eqn{b_a(t)}, \eqn{c_b(s)}, and \eqn{d_c(u)} are B-spline bases over the
predictor, response, and quantile domains. The intercept surface is represented
separately as
\deqn{
\alpha(s,u)
\approx
\sum_{b=1}^{L_s}
\sum_{c=1}^{L_u}
\xi_{bc}\widetilde c_b(s)\widetilde d_c(u).
}
The finite-dimensional linear predictor is therefore
\deqn{
\eta_i(s,u)
=
\{\widetilde d(u)\otimes \widetilde c(s)\}^{\top}\xi
+
\{d(u)\otimes c(s)\otimes r_i\}^{\top}\theta,
}
where
\deqn{
r_i =
\int_{\mathcal{T}} X_i(t)b(t)\,dt
}
is the predictor-basis score vector. The integral is computed numerically using
the rule specified by \code{integration}.

Estimation proceeds in three stages. First, a smooth pilot fit is obtained using
the empirical check loss together with directional P-spline smoothness penalties
and small ridge stabilizers. Second, an adaptive locally sparse estimator is
computed by adding a multiplicity-balanced overlapping group-Lasso penalty to
the coefficient tensor. This penalty acts on cuboidal blocks of B-spline
coefficients and is designed to identify local regions where the regression
kernel is exactly zero. Third, a constrained refit is performed on the selected
active coefficient set to reduce the shrinkage bias introduced by the group
penalty.

If \code{tuning = TRUE}, the function evaluates a finite collection of tuning
vectors and selects the one minimizing a QBIC-type criterion. The criterion is
computed after the full three-stage procedure and penalizes the effective model
dimension after support selection. If \code{tuning = FALSE}, the model is fitted
using \code{base_tuning} directly.

The fitted object stores the pilot, sparse, and refitted estimates. Prediction
and coefficient extraction use the final constrained-refit estimator
\code{fit_r}.
}

\value{
An object of class \code{"sfqfofr"}, represented as a list with the following
components:

\item{\code{fit_p}}{
  Pilot-stage fit. This is the smooth initial estimator used to construct
  adaptive group weights.
}

\item{\code{fit_s}}{
  Sparse-stage fit. This is the adaptive locally sparse estimator obtained
  before the constrained refit.
}

\item{\code{fit_r}}{
  Final constrained-refit estimator. This component is used by
  \code{\link{predict.sfqfofr}} and \code{\link{coef.sfqfofr}}.
}

\item{\code{tuning}}{
  Named list containing the tuning parameters used in the final fit.
}

\item{\code{tuning_table}}{
  If \code{tuning = TRUE}, a table containing the evaluated tuning candidates
  and their criterion values. If \code{tuning = FALSE}, this component is
  \code{NULL}.
}

\item{\code{design_train}}{
  Internal sparse design object constructed from the training data.
}

\item{\code{ginfo}}{
  Internal object containing the overlapping group information for the
  coefficient tensor.
}

\item{\code{M_grp}}{
  Internal group-membership matrix used by the overlapping group penalty.
}

\item{\code{bases}}{
  B-spline basis matrices used in the model fit.
}

\item{\code{grids}}{
  Predictor, response, and quantile grids used in the model fit.
}

\item{\code{integration}}{
  Numerical integration rule used to approximate functional inner products.
}
}
\author{
Ufuk Beyaztas
}
\note{
The input curves should be observed on the same grids used to construct the
basis matrices. If the data are observed on irregular or subject-specific grids,
they should first be registered or interpolated to common grids.

The default numerical integration rule is \code{"trapezoid"}, which is robust
for arbitrary grid lengths. Simpson integration can be requested explicitly
through \code{integration = "simpson"}, but it requires an odd number of
equally spaced grid points.

The function returns the final constrained-refit estimator as part of a classed
object. Use \code{\link{predict.sfqfofr}} for predicted conditional quantile
curves and \code{\link{coef.sfqfofr}} for estimated coefficient surfaces.
}
\examples{
\donttest{
## Construct grids
t_grid   <- seq(0, 1, length.out = 101)
s_grid   <- seq(0, 1, length.out = 101)
tau_grid <- c(0.10, 0.35, 0.50, 0.65, 0.80, 0.90)

## Build B-spline bases
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

## Generate a small example data set
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

## Predict conditional quantile response curves for test predictors
pred <- predict(fit, Xnew = dat$test$X)

## Dimensions: n_test by length(s_grid) by length(tau_grid)
dim(pred)

## Median quantile index
tau_id <- which.min(abs(tau_grid - 0.50))

## Select the first 10 test curves
id_plot <- 1:10

## True and estimated conditional median curves
Q_true_med <- dat$test$Qtrue[id_plot, , tau_id]
Q_hat_med  <- pred[id_plot, , tau_id]

## Plot true and estimated median curves
op <- par(
  mfrow = c(2, 5),
  mar = c(3.0, 3.0, 2.2, 0.8),
  oma = c(0.2, 0.2, 1.0, 0.2)
)

for (ii in seq_along(id_plot)) {

  yy <- range(c(Q_true_med[ii, ], Q_hat_med[ii, ]), na.rm = TRUE)

  plot(
    s_grid,
    Q_true_med[ii, ],
    type = "l",
    lwd = 2,
    ylim = yy,
    xlab = "s",
    ylab = "Median curve",
    main = paste0("Test curve ", id_plot[ii])
  )

  lines(
    s_grid,
    Q_hat_med[ii, ],
    lwd = 2,
    lty = 2
  )

  legend(
    "topright",
    legend = c("True", "Estimated"),
    lty = c(1, 2),
    lwd = c(2, 2),
    bty = "n",
    cex = 0.75
  )
}

mtext(
  expression("True and estimated conditional median curves, " * u == 0.50),
  side = 3,
  outer = TRUE,
  line = -0.2,
  font = 2
)

par(op)
}
}
