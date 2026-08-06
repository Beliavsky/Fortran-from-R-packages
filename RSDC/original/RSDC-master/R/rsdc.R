#' @name RSDC
#' @aliases RSDC-package
#' @title RSDC: Regime-Switching Dynamic Correlation Models
#'
#' @description
#' The \pkg{RSDC} package provides a comprehensive framework for modeling,
#' estimating, and forecasting correlation structures in multivariate time
#' series under regime-switching dynamics. It supports both fixed transition
#' probabilities and \emph{time-varying transition probabilities} (TVTP)
#' driven by exogenous variables.
#'
#' The methodology is particularly suited to empirical asset pricing and
#' portfolio management applications, enabling users to incorporate macroeconomic,
#' financial, or climate-related predictors into the regime dynamics. The package
#' integrates the full workflow — from model estimation to covariance matrix
#' reconstruction and portfolio optimization — in a single, reproducible pipeline.
#'
#' @section Main Features:
#' \itemize{
#'   \item \strong{Model estimation and filtering:}
#'     \code{\link{rsdc_hamilton}} (Hamilton filter),
#'     \code{\link{rsdc_likelihood}} (likelihood computation),
#'     \code{\link{rsdc_estimate}} (parameter estimation).
#'   \item \strong{Correlation and covariance forecasting:}
#'     \code{\link{rsdc_forecast}} (in-sample / 70-30 path),
#'     \code{\link{rsdc_forecast_ahead}} (multi-step-ahead),
#'     \code{\link{rsdc_corr_bands}} (parameter-uncertainty bands),
#'     \code{\link{rsdc_viterbi}} (most likely regime path).
#'   \item \strong{Portfolio construction:}
#'     \code{\link{rsdc_minvar}} (minimum-variance portfolios),
#'     \code{\link{rsdc_maxdiv}} (maximum-diversification portfolios).
#'   \item \strong{Simulation:}
#'     \code{\link{rsdc_simulate}} (simulate TVTP regime-switching series).
#'   \item \strong{Fitted-model methods:} \code{rsdc_estimate()} returns an
#'     object of class \code{"rsdc_fit"} with \code{print}, \code{summary},
#'     \code{coef}, \code{logLik} (so \code{AIC}/\code{BIC} work), \code{nobs},
#'     \code{vcov}, \code{confint}, \code{predict}, and \code{simulate} methods.
#' }
#'
#' @section Authors:
#' David Ardia, Benjamin Seguin and Roosevelt Ymele Nguemo. David Ardia is the
#' maintainer (also copyright holder and funder).
#'
#' @references
#' \insertRef{engle2002dynamic}{RSDC} \cr
#'
#' \insertRef{hamilton1989}{RSDC} \cr
#'
#' \insertRef{pelletier2006regime}{RSDC}
#'
#' @examples
#' \donttest{
#' # Quickstart: simulate two correlation regimes, fit, and inspect the fit
#' sim <- rsdc_simulate(n = 500, X = matrix(1, 500, 1),
#'                      beta = matrix(qlogis(0.9), 2, 1),
#'                      mu = matrix(0, 2, 2),
#'                      sigma = array(c(1, 0.1, 0.1, 1,
#'                                      1, 0.8, 0.8, 1), c(2, 2, 2)),
#'                      N = 2, seed = 2)
#' fit <- rsdc_estimate("noX", residuals = sim$observations, N = 2)
#' summary(fit)
#' c(AIC = AIC(fit), BIC = BIC(fit))
#' # Most likely regime path
#' table(rsdc_viterbi(fit))
#' # Getting-started vignette:  browseVignettes("RSDC")
#' # The companion article (methodology + a five-industry out-of-sample study)
#' # is reproducible from the packaged data; see citation("RSDC").
#' }
#'
#' @useDynLib RSDC, .registration = TRUE
#' @importFrom Rdpack reprompt
#' @importFrom Rcpp sourceCpp
#' @importFrom DEoptim DEoptim
#' @importFrom stats optim plogis cor sd
#' @importFrom stats nobs coef logLik vcov confint predict simulate
#' @importFrom utils combn
NULL

