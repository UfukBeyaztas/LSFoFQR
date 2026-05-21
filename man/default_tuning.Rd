\name{default_tuning}
\alias{default_tuning}
\title{Default Tuning Parameters for LSFoFQR}

\description{
Constructs a default tuning-parameter vector for the proposed locally sparse
simultaneous function-on-function quantile regression estimator. The returned
list contains the smoothing, sparsity, refitting, adaptivity, and support-screening
parameters used by \code{\link{sfqfofr}}.
}

\usage{
default_tuning(n_train)
}

\arguments{
  \item{n_train}{
    Positive integer giving the number of independent training curves. The value
    is used to mildly scale the group-sparsity parameter and to determine the
    default adaptive-weight exponent and active-block threshold.
  }
}

\details{
The proposed estimator is computed through a three-stage procedure. First, a
smooth pilot fit is obtained using only smoothness and ridge regularization.
Second, an adaptive locally sparse estimator is computed by adding a
multiplicity-balanced overlapping group-Lasso penalty to the coefficient tensor.
Third, a constrained refit is performed on the active coefficient set selected
from the sparse fit.

The tuning vector returned by \code{default_tuning()} contains the numerical
parameters controlling these three stages. The pilot-stage smoothing parameters
are \code{lam_init_b} for the coefficient kernel and \code{lam_init_a} for the
intercept surface. The sparse-stage smoothing parameters are \code{lam_sp_b}
and \code{lam_sp_a}, while \code{lam_sp_g} controls the adaptive overlapping
group-Lasso penalty. The refit-stage smoothing parameters are \code{lam_rf_b}
and \code{lam_rf_a}. The support-screening step is controlled by
\code{active_tol}, \code{active_method}, and \code{min_gap}. The adaptive-weight
exponent is given by \code{eta_adapt}.

The simultaneous function-on-function quantile regression model is
\deqn{
Q_{Y(s)\mid X}(u \mid X)
=
\alpha(s,u)
+
\int_{\mathcal{T}} X(t)\beta(t,s,u)\,dt,
}
where \eqn{\beta(t,s,u)} is represented by a tensor-product B-spline expansion
over the predictor, response, and quantile domains. The local-sparsity penalty
acts on overlapping cuboidal blocks of the coefficient tensor. Therefore, the
tuning parameters returned by this function determine the balance between
smoothness, sparsity, and post-selection refitting in the estimated kernel
\eqn{\widehat{\beta}(t,s,u)}.

The default values are intended to provide stable initial values for practical
use and for construction of tuning grids. They may not be optimal for every data
set. Users can modify the returned list directly or use it as the center of a
candidate grid through \code{\link{make_tuning_grid}}.
}

\value{
A named list with the following components:

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
  Sparse-stage adaptive overlapping group-Lasso parameter. Larger values produce
  stronger local sparsity in the estimated coefficient kernel. The default value
  is scaled by \eqn{\sqrt{n_{\mathrm{train}}/100}}.
}

\item{\code{lam_rf_b}}{
  Refit-stage smoothing parameter for the coefficient kernel after support
  selection.
}

\item{\code{lam_rf_a}}{
  Refit-stage smoothing parameter for the intercept surface after support
  selection.
}

\item{\code{active_tol}}{
  Threshold used in the hybrid active-block screening rule. Blocks with
  sufficiently large hybrid support scores are retained for the constrained
  refit. The default is bounded between \code{0.03} and \code{0.45}.
}

\item{\code{eta_adapt}}{
  Adaptive-weight exponent used in the group penalty. Larger values place more
  emphasis on the pilot fit when constructing adaptive group weights.
}

\item{\code{active_method}}{
  Character string specifying the active-block selection rule. The default
  \code{"fixed_or_gap"} combines a fixed threshold with a gap-based fallback.
}

\item{\code{min_gap}}{
  Minimum gap required for the gap-based active-block screening rule.
}
}
\author{
Ufuk Beyaztas
}

\note{
The values returned by this function are default starting values, not universal
optimal tuning parameters. For formal analysis, the tuning parameters should be
selected by the QBIC-type criterion implemented in the fitting procedure or by a
user-defined tuning grid.
}
\examples{
## Default tuning for a sample of 100 functional observations
tun100 <- default_tuning(n_train = 100)
tun100

## Default tuning for a larger sample
tun250 <- default_tuning(n_train = 250)
tun250

## Modify one component manually
tun100$lam_sp_g <- 2 * tun100$lam_sp_g

## Use the default vector as the center of a tuning grid
grid <- make_tuning_grid(n_train = 100, base_tuning = tun100)
length(grid)
}
