#' @title AIC for Tweedie Glms
#' @name tweedie_AIC
#' @description Evaluates the \acronym{aic} for a fitted Tweedie \acronym{glm}.
#' The Tweedie family of distributions belong to the class of exponential dispersion models (\acronym{edm}s), 
#' famous for their role in generalized linear models. 
#' The Tweedie distributions are the \acronym{edm}s with a variance of the form 
#' \eqn{\mbox{var}[Y] = \phi\mu^p}{var[Y] = phi*mu^p} where \eqn{p \ge 1}{p >= 1}.
#' \emph{This function only evaluates for \eqn{p \ge 1}{p >= 1}.}
#'
#' @usage tweedie_AIC(glm.obj, dispersion = NULL, k = 2, verbose = TRUE)
#' 
#' @details
#' The \acronym{aic} is computed by evaluating the density function.

#' @note
#' Evaluating the likelihood can be time consuming, so the function may take some time for large data sets.
#'
#' @param glm.obj a fitted \code{glm} object, fitted using the \code{tweedie} family.
#' @param dispersion the dispersion parameter, usually extracted from \code{glm.obj}; however, occasionally  a specified value of the dispersion may be needed.
#' @param k the \acronym{aic} penalty; \code{k = 2} (the default) produces the AIC.
#' @param verbose logical; if \code{TRUE}, display details of the internal process. The default is \code{FALSE}.
#'
#' @return The value of the computed \acronym{aic}.
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
#' tweedie_AIC(fit)
#'
#' @export
tweedie_AIC <- function( glm.obj, dispersion = NULL, k = 2, verbose = TRUE){ 
  # New  dispersion  input for (e.g.) Poisson case, added October 2017
  
  wt <- glm.obj$prior.weights
  n <- length(glm.obj$residuals)
  edf <- glm.obj$rank  # As used in logLik.glm()
  
  mu <- stats::fitted( glm.obj )
  y  <- glm.obj$y
  p <- get("p", envir = environment(glm.obj$family$variance))
  
  if ( is.null(dispersion)) {  # New section
    if (p == 1 & verbose) message("*** Tweedie index power = 1: Consider using  dispersion=1  in call to  tweedie_AIC().\n")
    dev <- deviance(glm.obj)
    disp <- dev / sum(wt)  # In line with Gamma()$aic
    edf <- edf + 1  # ADD one as we are estimating phi too
  } else {
    disp <- dispersion
  }
  
  den <- dtweedie( y = y, 
                   mu = mu, 
                   phi = disp, 
                   power = p)
  AIC <- -2 * sum( log(den) * wt) 
  
  return( AIC + k * (edf) )
  
}


#' @rdname tweedie_AIC
#' @export
AICtweedie <- function( glm.obj, dispersion = NULL, k = 2, verbose = TRUE){ 
  lifecycle::deprecate_warn(when = "3.0.5", 
                            what = "AICtweedie()", 
                            with = "tweedie_AIC()")
  tweedie_AIC(glm.obj = glm.obj, 
              dispersion = dispersion,
              k = k,
              verbose = verbose)
}


