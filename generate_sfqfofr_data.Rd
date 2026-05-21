\name{generate_sfqfofr_data}
\alias{generate_sfqfofr_data}
\title{Generate Simulated SFQFoFR Data}

\description{
Generates simulated functional predictor and functional response data for the
simultaneous function-on-function quantile regression model. The data-generating
mechanism follows a location-scale functional model with a known
quantile-varying coefficient surface, allowing the fitted coefficient kernel and
conditional quantile surfaces to be evaluated against the truth.
}

\usage{
generate_sfqfofr_data(
  n_train = 100,
  n_test = 1000,
  t_grid,
  s_grid,
  tau_grid,
  bases,
  error_type = c("gaussian", "t3", "cauchy"),
  support_design = c("nonsparse", "nonsparse_noncrossing", "sparse"),
  integration = c("trapezoid", "simpson", "auto")
)
}

\arguments{
  \item{n_train}{
    Positive integer giving the number of training functional observations.
  }

  \item{n_test}{
    Positive integer giving the number of test functional observations.
  }

  \item{t_grid}{
    Numeric vector giving the grid over the predictor domain \eqn{\mathcal{T}}.
    The grid must be strictly increasing.
  }

  \item{s_grid}{
    Numeric vector giving the grid over the response domain \eqn{\mathcal{S}}.
    The grid must be strictly increasing.
  }

  \item{tau_grid}{
    Numeric vector of quantile levels. All elements should lie strictly between
    0 and 1.
  }

  \item{bases}{
    A list of B-spline basis matrices, typically returned by
    \code{\link{make_sfqfofr_basis}}. The list must contain the components
    \code{Bt}, \code{Bs}, \code{Bu}, \code{Bsa}, and \code{Bua}.
  }

  \item{error_type}{
    Character string specifying the marginal innovation distribution. Available
    options are \code{"gaussian"}, \code{"t3"}, and \code{"cauchy"}.
    The default is \code{"gaussian"}.
  }

  \item{support_design}{
    Character string specifying the coefficient-surface design. Available
    options are \code{"nonsparse"}, \code{"nonsparse_noncrossing"}, and
    \code{"sparse"}.
  }

  \item{integration}{
    Numerical integration rule used to approximate functional inner products.
    Available options are \code{"trapezoid"}, \code{"simpson"}, and
    \code{"auto"}. The default is \code{"trapezoid"}, which is recommended for
    general grids because it works for arbitrary grid lengths.
  }
}

\details{
The function generates data from a functional location-scale model whose
conditional quantile function has the simultaneous function-on-function
quantile regression form
\deqn{
Q_{Y(s)\mid X}(u \mid X)
=
\alpha_0(s,u)
+
\int_{\mathcal{T}} X(t)\beta_0(t,s,u)\,dt.
}
The functional predictor is generated from a finite Karhunen--Loeve-type
expansion,
\deqn{
X_i(t)
=
\sum_{k=1}^{6} \xi_{ik}\phi_k(t) + \delta_i(t),
}
where the basis functions are trigonometric functions and the scores have
decreasing standard deviations. The generated predictor trajectories are then
centered pointwise.

The functional response is generated from the location-scale model
\deqn{
Y_i(s)
=
\mu_i(s) + \sigma_i(s)\varepsilon_i(s),
}
where
\deqn{
\mu_i(s)
=
a_\mu(s)
+
\int_{\mathcal{T}} X_i(t)\Theta_\mu(t,s)\,dt,
}
and
\deqn{
\sigma_i(s)
=
a_\sigma(s)
+
\int_{\mathcal{T}} X_i(t)\Theta_\sigma(t,s)\,dt.
}
The innovation process \eqn{\varepsilon_i(s)} has exponential correlation over
the response domain. The marginal innovation distribution is selected by
\code{error_type}. The Gaussian setting provides a light-tailed case, the
Student-\eqn{t_3} setting provides a moderately heavy-tailed case, and the
Cauchy setting provides an extreme heavy-tailed stress test.

Because the innovation distribution is symmetric, the conditional quantile
surface is available in closed form:
\deqn{
Q_{Y_i(s)\mid X_i}(u\mid X_i)
=
a_\mu(s) + q_u a_\sigma(s)
+
\int_{\mathcal{T}} X_i(t)
\left\{
\Theta_\mu(t,s) + q_u\Theta_\sigma(t,s)
\right\}
dt,
}
where \eqn{q_u} is the \eqn{u}th quantile of the selected innovation
distribution. Consequently,
\deqn{
\alpha_0(s,u)=a_\mu(s)+q_u a_\sigma(s),
\qquad
\beta_0(t,s,u)=\Theta_\mu(t,s)+q_u\Theta_\sigma(t,s).
}
The returned object contains both the observed training and test samples and the
corresponding true conditional quantile surfaces.

Three coefficient-support designs are available. The \code{"nonsparse"} design
uses a dense smooth coefficient surface that changes sign. The
\code{"nonsparse_noncrossing"} design uses a dense smooth coefficient surface
that is constrained to remain nonzero over the considered quantile range. The
\code{"sparse"} design generates a locally sparse coefficient surface with
separated active regions and extended inactive regions. The sparse design
requires at least 11 basis functions in both the predictor and response
directions.
}

\value{
A list with the following components:

\item{train}{
  A list containing the training data:
  \itemize{
    \item \code{X}: an \code{n_train} by \code{length(t_grid)} matrix of
    functional predictor observations;
    \item \code{Y}: an \code{n_train} by \code{length(s_grid)} matrix of
    functional response observations;
    \item \code{Qtrue}: an array of dimension
    \code{n_train} by \code{length(s_grid)} by \code{length(tau_grid)}
    containing the true conditional quantile response surfaces for the training
    observations.
  }
}

\item{test}{
  A list containing the test data:
  \itemize{
    \item \code{X}: an \code{n_test} by \code{length(t_grid)} matrix of
    functional predictor observations;
    \item \code{Y}: an \code{n_test} by \code{length(s_grid)} matrix of
    functional response observations;
    \item \code{Qtrue}: an array of dimension
    \code{n_test} by \code{length(s_grid)} by \code{length(tau_grid)}
    containing the true conditional quantile response surfaces for the test
    observations.
  }
}

\item{truth}{
  A list containing the true data-generating objects:
  \itemize{
    \item \code{Theta_mean}: coefficient matrix for the location component;
    \item \code{Theta_scale}: coefficient matrix for the scale component;
    \item \code{theta_true_arr}: tensor of true B-spline coefficients;
    \item \code{alpha_coef}: matrix of true intercept-surface coefficients;
    \item \code{beta_arr}: array containing the true coefficient surfaces
    \eqn{\beta_0(t,s,u)} on the observed grids;
    \item \code{alpha_arr}: matrix containing the true intercept surface
    \eqn{\alpha_0(s,u)};
    \item \code{mean_curve}: generated conditional location curves;
    \item \code{scale_curve}: generated conditional scale curves;
    \item \code{q_tau}: distributional quantiles corresponding to
    \code{tau_grid};
    \item \code{support_design}: support design used in the simulation;
    \item \code{error_type}: innovation distribution used in the simulation.
  }
}
}
\author{
Ufuk Beyaztas
}
\examples{
\donttest{

t_grid <- seq(0, 1, length.out = 101)
s_grid <- seq(0, 1, length.out = 101)
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

dim(dat$train$X)
dim(dat$train$Y)
dim(dat$train$Qtrue)
dim(dat$truth$beta_arr)
}
}
