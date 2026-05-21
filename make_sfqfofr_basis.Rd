\name{make_sfqfofr_basis}
\alias{make_sfqfofr_basis}
\title{Construct B-Spline Bases for SFQFoFR}

\description{
Constructs the B-spline basis matrices used to fit the locally sparse
simultaneous function-on-function quantile regression model. The function
creates separate bases for the predictor domain, response domain, quantile
domain, and the intercept surface.
}

\usage{
make_sfqfofr_basis(
  t_grid,
  s_grid,
  tau_grid,
  kt = 11,
  ks = 11,
  ku = NULL,
  kas = 5,
  kau = NULL,
  degree_t = 2L,
  degree_s = 2,
  degree_u = 2
)
}

\arguments{
  \item{t_grid}{
    Numeric vector giving the grid over the predictor domain
    \eqn{\mathcal{T}}. The grid should be strictly increasing.
  }

  \item{s_grid}{
    Numeric vector giving the grid over the response domain
    \eqn{\mathcal{S}}. The grid should be strictly increasing.
  }

  \item{tau_grid}{
    Numeric vector of quantile levels over the quantile domain
    \eqn{\mathcal{U}}. All entries should lie strictly between 0 and 1.
  }

  \item{kt}{
    Number of B-spline basis functions used in the predictor direction
    \eqn{t}. The default is \code{11}.
  }

  \item{ks}{
    Number of B-spline basis functions used in the response direction
    \eqn{s}. The default is \code{11}.
  }

  \item{ku}{
    Number of B-spline basis functions used in the quantile direction
    \eqn{u} for the coefficient kernel. If \code{NULL}, it is set equal to
    \code{length(tau_grid)}.
  }

  \item{kas}{
    Number of B-spline basis functions used in the response direction
    for the intercept surface \eqn{\alpha(s,u)}. The default is \code{5}.
  }

  \item{kau}{
    Number of B-spline basis functions used in the quantile direction
    for the intercept surface \eqn{\alpha(s,u)}. If \code{NULL}, it is set
    equal to \code{length(tau_grid)}.
  }

  \item{degree_t}{
    Polynomial degree of the B-spline basis in the predictor direction.
    The default is \code{2}, corresponding to quadratic B-splines.
  }

  \item{degree_s}{
    Polynomial degree of the B-spline basis in the response direction.
    The default is \code{2}, corresponding to quadratic B-splines.
  }

  \item{degree_u}{
    Polynomial degree of the B-spline basis in the quantile direction.
    The default is \code{2}, corresponding to quadratic B-splines.
  }
}

\details{
The simultaneous function-on-function quantile regression model has the form
\deqn{
Q_{Y(s)\mid X}(u \mid X)
=
\alpha(s,u)
+
\int_{\mathcal{T}} X(t)\beta(t,s,u)\,dt,
\qquad (s,u)\in\mathcal{S}\times\mathcal{U}.
}
The coefficient kernel \eqn{\beta(t,s,u)} is represented by a tensor-product
B-spline expansion,
\deqn{
\beta(t,s,u)
\approx
\sum_{a=1}^{K_t}
\sum_{b=1}^{K_s}
\sum_{c=1}^{K_u}
\theta_{abc} b_a(t)c_b(s)d_c(u),
}
where \eqn{b_a(t)}, \eqn{c_b(s)}, and \eqn{d_c(u)} are B-spline bases over
the predictor, response, and quantile domains, respectively.

The intercept surface is represented separately as
\deqn{
\alpha(s,u)
\approx
\sum_{b=1}^{L_s}
\sum_{c=1}^{L_u}
\xi_{bc}\widetilde c_b(s)\widetilde d_c(u).
}
This function constructs the basis matrices required for these two
representations. The coefficient-kernel basis consists of \code{Bt}, \code{Bs},
and \code{Bu}. The intercept-surface basis consists of \code{Bsa} and
\code{Bua}.

The tensor-product representation is used internally by \code{\link{sfqfofr}}
to build the finite-dimensional linear predictor,
\deqn{
\eta_i(s,u)
=
\{\widetilde d(u)\otimes \widetilde c(s)\}^\top \xi
+
\{d(u)\otimes c(s)\otimes r_i\}^\top \theta,
}
where
\deqn{
r_i =
\int_{\mathcal{T}} X_i(t)b(t)\,dt
}
is the predictor-basis score vector. The local support of the B-spline basis is
also central to the local-sparsity mechanism: when a block of tensor-product
basis coefficients is set to zero, the fitted coefficient kernel is exactly zero
on the corresponding open tensor-product knot cells.

For the default quadratic B-splines, each basis function has compact local
support over a small number of adjacent intervals. This property allows the
overlapping group-Lasso penalty in \code{\link{sfqfofr}} to induce interpretable
locally inactive regions in the estimated coefficient surface.
}

\value{
A named list with the following components:

\item{\code{Bt}}{
  Matrix of B-spline basis evaluations over \code{t_grid}. It has dimension
  \code{length(t_grid)} by \code{kt}. This basis is used for the predictor
  direction of the coefficient kernel.
}

\item{\code{Bs}}{
  Matrix of B-spline basis evaluations over \code{s_grid}. It has dimension
  \code{length(s_grid)} by \code{ks}. This basis is used for the response
  direction of the coefficient kernel.
}

\item{\code{Bu}}{
  Matrix of B-spline basis evaluations over \code{tau_grid}. It has dimension
  \code{length(tau_grid)} by \code{ku}. This basis is used for the quantile
  direction of the coefficient kernel.
}

\item{\code{Bsa}}{
  Matrix of B-spline basis evaluations over \code{s_grid}. It has dimension
  \code{length(s_grid)} by \code{kas}. This basis is used for the response
  direction of the intercept surface.
}

\item{\code{Bua}}{
  Matrix of B-spline basis evaluations over \code{tau_grid}. It has dimension
  \code{length(tau_grid)} by \code{kau}. This basis is used for the quantile
  direction of the intercept surface.
}
}
\author{
Ufuk Beyaztas
}

\note{
The default polynomial degree is \code{2}, corresponding to quadratic B-splines.
The locally sparse estimator uses overlapping cuboidal blocks of basis
coefficients. For stable use of the default local-sparsity structure, the basis
dimensions \code{kt}, \code{ks}, and \code{ku} should generally be at least
\code{3}.

When \code{ku = NULL} or \code{kau = NULL}, the corresponding number of
quantile-basis functions is set equal to \code{length(tau_grid)}. This is
convenient when the model is fitted on a modest fixed grid of quantile levels.
}
\examples{
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

## Inspect basis dimensions
dim(bases$Bt)
dim(bases$Bs)
dim(bases$Bu)
dim(bases$Bsa)
dim(bases$Bua)

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


## Use the basis object in a small simulated example
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
