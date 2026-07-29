#' bayesianOU: Bayesian Nonlinear Ornstein-Uhlenbeck Models
#'
#' Fits Bayesian nonlinear Ornstein-Uhlenbeck models with cubic drift,
#' stochastic volatility (SV), and Student-t innovations. Sector-specific
#' parameters use a non-centered hierarchical (partial-pooling) parameterization,
#' and sampling runs via 'Stan' (with within-chain parallelism through
#' \code{reduce_sum}).
#'
#' @section Main Functions:
#' \itemize{
#'    \item \code{\link{fit_ou_nonlinear_tmg}}: Fit the main OU model
#'    \item \code{\link{extract_posterior_summary}}: Extract posterior summaries
#'    \item \code{\link{validate_ou_fit}}: Validate model fit (MCMC + LOO + OOS)
#'    \item \code{\link{kappa_stability_evidence}}: Dynamic mean-reversion evidence
#'    \item \code{\link{compare_models_loo}}: Compare models via PSIS-LOO
#' }
#'
#' @section Model Specification:
#' The estimated one-step (Euler, dt = 1) mean equation is
#' \deqn{\Delta Y_{t,s} = \kappa_s(\theta_s - Y_{t-1,s} + a_3 (Y_{t-1,s}-\theta_s)^3)
#'   + (\beta_{0,s} + \beta_1 zTMG_t) X_{t-1,s} + \gamma\, COM_{t-1,s} + \epsilon_{t,s}}
#' with \eqn{\epsilon \sim \mathrm{Student\text{-}t}(\nu, 0, \sigma_{t,s})} and
#' \eqn{\log\sigma^2_{t,s}} a stationary AR(1). See the README for the full
#' methodology (train/test design, priors, and the \eqn{a_3<0} assumption).
#'
#' @docType package
#' @name bayesianOU-package
#' @aliases bayesianOU
#' @keywords internal
"_PACKAGE"

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "bayesianOU: Bayesian Nonlinear OU Models\n",
    "For 'Stan' backend, ensure 'cmdstanr' or 'rstan' is installed."
  )
}