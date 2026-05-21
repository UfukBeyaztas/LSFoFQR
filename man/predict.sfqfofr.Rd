\name{predict.sfqfofr}
\alias{predict.sfqfofr}
\alias{predict}
\title{Predict Conditional Quantile Curves from an \code{sfqfofr} Fit}

\description{
Computes predicted conditional quantile response curves for new functional
predictors using a fitted locally sparse simultaneous function-on-function
quantile regression model.
}

\usage{
\method{predict}{sfqfofr}(object, Xnew, ...)
}

\arguments{
  \item{object}{
    A fitted model object of class \code{"sfqfofr"}, typically returned by
    \code{\link{sfqfofr}}.
  }

  \item{Xnew}{
    Numeric matrix containing new functional predictor observations. Each row
    represents one new functional predictor curve, and each column corresponds
    to one point on the predictor grid \code{object$grids$t}. Thus,
    \code{ncol(Xnew)} must equal \code{length(object$grids$t)}.
  }

  \item{...}{
    Additional arguments passed to the S3 method. Currently unused.
  }
}

\details{
The fitted model has the simultaneous function-on-function quantile regression
form
\deqn{
Q_{Y(s)\mid X}(u \mid X)
=
\alpha(s,u)
+
\int_{\mathcal{T}} X(t)\beta(t,s,u)\,dt,
\qquad (s,u)\in\mathcal{S}\times\mathcal{U}.
}
For a new functional predictor \eqn{X^\star(t)}, this function computes
\deqn{
\widehat Q_{Y^\star(s)\mid X^\star}(u \mid X^\star)
=
\widehat\alpha(s,u)
+
\int_{\mathcal{T}} X^\star(t)\widehat\beta(t,s,u)\,dt
}
on the response grid and quantile grid stored in the fitted object.

Internally, the functional inner product is approximated by first computing the
predictor-basis score vector
\deqn{
\widehat r^\star
=
\sum_j \delta_j X^\star(t_j)b(t_j),
}
where \eqn{\delta_j} are numerical integration weights and \eqn{b(t)} is the
B-spline basis over the predictor domain. The integration rule is inherited
from the fitted object through \code{object$integration}. The default fitting
rule is trapezoidal integration, which is suitable for general grids.

The coefficient kernel used in prediction is the final constrained-refit
estimator stored in \code{object$fit_r}. Therefore, predicted quantile curves
are based on the locally sparse, bias-reduced final estimator rather than on
the pilot or sparse intermediate fits.

The prediction is simultaneous over all quantile levels in
\code{object$grids$tau}. The output is an array whose third dimension indexes
the quantile levels.
}

\value{
A three-dimensional numeric array with dimensions
\deqn{
n_{\mathrm{new}}
\times
\texttt{length(object\$grids\$s)}
\times
\texttt{length(object\$grids\$tau)}.
}
The first dimension indexes new observations, the second dimension indexes the
response-domain grid, and the third dimension indexes the quantile levels.

For example, \code{out[i, , r]} is the predicted conditional quantile response
curve for the \code{i}th new predictor at the quantile level
\code{object$grids$tau[r]}.
}
\author{
Ufuk Beyaztas
}
\note{
The new predictor matrix \code{Xnew} must be observed on the same predictor
grid used to fit the model. If new curves are observed on a different grid, they
should first be interpolated or registered to \code{object$grids$t}.

This function returns predicted conditional quantile curves. To extract the
estimated coefficient kernel \eqn{\widehat\beta(t,s,u)}, use
\code{\link{coef.sfqfofr}}.
}
\examples{
\donttest{
## Construct grids
t_grid   <- seq(0, 1, length.out = 101)
s_grid   <- seq(0, 1, length.out = 101)
tau_grid <- c(0.10, 0.35, 0.50, 0.65, 0.80, 0.90)

## Build basis matrices
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

## Generate a small data set
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
