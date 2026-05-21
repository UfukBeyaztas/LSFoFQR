\name{make_tuning_grid}
\alias{make_tuning_grid}
\title{Construct a Tuning Grid for LSFoFQR}

\description{
Constructs a finite candidate grid of tuning-parameter vectors for the proposed
locally sparse simultaneous function-on-function quantile regression estimator.
The grid is formed by multiplying a baseline tuning vector by user-specified
smoothness and sparsity factors and by varying the active-block screening
parameters.
}

\usage{
make_tuning_grid(
  n_train,
  base_tuning = default_tuning(n_train),
  smooth_mult = c(0.5, 1.0, 2.0),
  group_mult = c(0.5, 1.0, 2.0),
  active_tol_grid = NULL,
  eta_grid = NULL,
  active_method = "fixed_or_gap",
  min_gap = 0.02
)
}

\arguments{
  \item{n_train}{
    Positive integer giving the number of independent training curves. This
    argument is used to construct \code{base_tuning} through
    \code{\link{default_tuning}} when \code{base_tuning} is not supplied.
  }

  \item{base_tuning}{
    Named list giving the baseline tuning vector around which the grid is
    constructed. By default, this is \code{default_tuning(n_train)}.
  }

  \item{smooth_mult}{
    Numeric vector of multiplicative factors applied to all smoothness
    parameters. These factors are applied to \code{lam_init_b},
    \code{lam_init_a}, \code{lam_sp_b}, \code{lam_sp_a}, \code{lam_rf_b}, and
    \code{lam_rf_a}.
  }

  \item{group_mult}{
    Numeric vector of multiplicative factors applied to the adaptive
    overlapping group-Lasso parameter \code{lam_sp_g}.
  }

  \item{active_tol_grid}{
    Optional numeric vector of active-block threshold values. If \code{NULL},
    the grid is constructed automatically around \code{base_tuning$active_tol}
    using the values \code{base_tuning$active_tol - 0.08},
    \code{base_tuning$active_tol}, and \code{base_tuning$active_tol + 0.08},
    truncated to lie between \code{0.03} and \code{0.55}.
  }

  \item{eta_grid}{
    Optional numeric vector of adaptive-weight exponents. If \code{NULL}, the
    value \code{base_tuning$eta_adapt} is used.
  }

  \item{active_method}{
    Character string specifying the support-screening rule used in the
    constrained refit. The default is \code{"fixed_or_gap"}.
  }

  \item{min_gap}{
    Positive numeric value giving the minimum gap required by the gap-based
    active-block screening rule. The default is \code{0.02}.
  }
}

\details{
The proposed estimator is computed through three stages: a smooth pilot fit, an
adaptive locally sparse fit, and a constrained refit on the selected active
coefficient set. The tuning parameters control smoothness, local sparsity,
adaptive weighting, and support screening across these stages.

Let
\deqn{
\zeta =
(\lambda_{\alpha}^{(0)},\lambda_{\beta}^{(0)},
 \lambda_{\alpha},\lambda_{\beta},\tau_{\beta},
 \lambda_{\alpha}^{\star},\lambda_{\beta}^{\star},
 \kappa,\eta)
}
denote a generic tuning vector. Here
\eqn{\lambda_{\alpha}^{(0)}} and \eqn{\lambda_{\beta}^{(0)}} control the smooth
pilot fit, \eqn{\lambda_{\alpha}}, \eqn{\lambda_{\beta}}, and
\eqn{\tau_{\beta}} control the adaptive sparse fit,
\eqn{\lambda_{\alpha}^{\star}} and \eqn{\lambda_{\beta}^{\star}} control the
constrained refit, \eqn{\kappa} is the active-block threshold, and \eqn{\eta}
is the adaptive-weight exponent.

This function constructs a finite candidate set \eqn{\mathcal{Z}} by scaling a
baseline tuning vector. The smoothness multiplier \code{smooth_mult} is applied
simultaneously to the pilot-stage, sparse-stage, and refit-stage smoothing
parameters. The group multiplier \code{group_mult} is applied only to the
adaptive overlapping group-Lasso parameter \code{lam_sp_g}. The threshold and
adaptive exponent are varied through \code{active_tol_grid} and
\code{eta_grid}.

The resulting list is designed for use in \code{\link{sfqfofr}} when
\code{tuning = TRUE}. For each candidate tuning vector, the full three-stage
estimation procedure is carried out, and the selected tuning vector is the one
minimizing the QBIC-type criterion implemented internally by the fitting
procedure.
}

\value{
A list of tuning-parameter vectors. Each element of the returned list is itself
a named list with the following components:

\item{\code{lam_init_b}}{
  Pilot-stage smoothing parameter for the coefficient kernel
  \eqn{\beta(t,s,u)}.
}

\item{\code{lam_init_a}}{
  Pilot-stage smoothing parameter for the intercept surface
  \eqn{\alpha(s,u)}.
}

\item{\code{lam_sp_b}}{
  Sparse-stage smoothing parameter for the coefficient kernel.
}

\item{\code{lam_sp_a}}{
  Sparse-stage smoothing parameter for the intercept surface.
}

\item{\code{lam_sp_g}}{
  Sparse-stage adaptive overlapping group-Lasso parameter controlling the
  strength of local sparsity.
}

\item{\code{lam_rf_b}}{
  Refit-stage smoothing parameter for the coefficient kernel.
}

\item{\code{lam_rf_a}}{
  Refit-stage smoothing parameter for the intercept surface.
}

\item{\code{active_tol}}{
  Threshold used in the hybrid active-block screening score.
}

\item{\code{eta_adapt}}{
  Adaptive-weight exponent used to construct the group-specific adaptive
  weights from the pilot fit.
}

\item{\code{active_method}}{
  Character string specifying the active-block screening rule.
}

\item{\code{min_gap}}{
  Minimum gap used by the gap-based part of the active-block screening rule.
}
}
\author{
Ufuk Beyaztas
}
\note{
If \code{base_tuning} is supplied directly, \code{n_train} is not used inside
the function body. However, when \code{base_tuning} is omitted, \code{n_train}
must be provided so that \code{\link{default_tuning}} can construct the baseline
tuning vector.
}
\examples{
## Construct a default tuning grid for n = 100
grid1 <- make_tuning_grid(n_train = 100)
length(grid1)

## Inspect the first candidate tuning vector
grid1[[1]]

## Use a smaller grid for faster tuning
grid2 <- make_tuning_grid(
  n_train = 100,
  smooth_mult = c(1),
  group_mult = c(0.5, 1, 2),
  active_tol_grid = c(0.15, 0.25),
  eta_grid = c(1.25)
)
length(grid2)

## Start from a user-modified baseline tuning vector
base <- default_tuning(n_train = 100)
base$lam_sp_g <- 2 * base$lam_sp_g

grid3 <- make_tuning_grid(
  n_train = 100,
  base_tuning = base,
  smooth_mult = c(0.5, 1, 2),
  group_mult = c(1)
)
length(grid3)
}
