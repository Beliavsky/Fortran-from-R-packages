#' @title Log-likelihood for Tweedie distributions
#' @name logLiktweedie
#' @description Evaluates the log-likelihood for a fitted Tweedie \acronym{glm}.
#'
#' @usage logLiktweedie(glm.obj, dispersion = NULL)
#' 
#' @details
#' The log-Likelihood is computed by evaluating the density function.

#' @note
#' Evaluating the likelihood can be time consuming, so the function may take some time for large data sets.
#'
#' @param glm.obj a fitted \code{glm} object, fitted using the \code{tweedie} family.
#' @param dispersion the dispersion parameter, usually extracted from \code{glm.obj}; however, occasionally  a specified value of the dispersion may be needed.
#'
#' @return The value of the computed log-likelihood.
#'
#' @seealso \code{\link{dtweedie}}
#'
#' @references
#' Dunn, P. K. and Smyth, G. K. (2008).
#' Evaluation of Tweedie exponential dispersion model densities by Fourier inversion.
#' \emph{Statistics and Computing}, 
#' \bold{18}, 73--86.
#' \doi{10.1007/s11222-007-9039-6}
#' 
#' Dunn, Peter K and Smyth, Gordon K (2005).
#' Series evaluation of Tweedie exponential dispersion model densities
#' \emph{Statistics and Computing},
#' \bold{15}(4). 267--280.
#' \doi{10.1007/s11222-005-4070-y}
#' 
#' Jorgensen, B. (1997).
#' \emph{Theory of Dispersion Models}.
#' Chapman and Hall, London.
#' 
#' Sakamoto, Y., Ishiguro, M., and Kitagawa G. (1986). 
#' \emph{Akaike Information Criterion Statistics}. 
#' D. Reidel Publishing Company.
#' 
#' @importFrom stats glm
#' @importFrom statmod tweedie
#'
#' @examples
#' # Fit a Tweedie density using  tweedie  family function from  statmod
#' pretend <- data.frame( y = stats::rgamma(20, shape = 1, rate = 1) )
#' fit <- glm(y ~ 1, data = pretend, 
#'            family = statmod::tweedie(link.power = 0, var.power = 2.1))
#'
#' # Compute the AIC
#' logLiktweedie(fit)
#'
#' @export
logLiktweedie <- function(glm.obj, dispersion = NULL) {
  # Computes the log-likelihood for a Tweedie glm.
  
  # Peter Dunn
  # 19 October 2017
  
  p <- get("p", envir = environment(glm.obj$family$variance))
  if (p == 1) message("*** Tweedie index power = 1: Consider using  dispersion = 1  in call to  logLiktweedie().\n")
  
  tweedie_AIC(glm.obj, 
              dispersion = dispersion, 
              k = 0, 
              verbose = FALSE) / (-2)
}
