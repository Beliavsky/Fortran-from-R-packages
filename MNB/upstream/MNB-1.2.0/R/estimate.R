l.MNB <- function(star, formula, dataSet) {

  Y <- stats::model.response(data = stats::model.frame(formula,dataSet))
  X <- stats::model.matrix(formula, dataSet)
  off <- stats::model.extract(stats::model.frame(formula,dataSet),"offset")

  dataSet.ind <- split(dataSet, f = dataSet$ind)
  n <- length(dataSet.ind)
  p <- dim(X)[2]
  mi <- dim(dataSet.ind[[1]])[1]
  N <- n * mi
  YI <- apply(matrix(Y, mi, n), 2, sum)

  phi <- star[1]
  beta <- star[2:(p + 1)]

  if(methods::is(off,"NULL")){
    eta <- X %*% beta
    mu.i <- apply(matrix(exp(eta), mi, n), 2, sum)
    logver <- sum(lgamma(phi + YI)) - n * lgamma(phi) -
      sum(lfactorial(Y)) + n * phi * (log(phi)) +
      t(Y) %*% eta - sum((phi + YI) * log(phi +
                                            mu.i))
    return(logver)}
  else{
    eta <- X %*% beta + off
    mu.i <- apply(matrix(exp(eta), mi, n), 2, sum)
    logver <- sum(lgamma(phi + YI)) - n * lgamma(phi) -
      sum(lfactorial(Y)) + n * phi * (log(phi)) +
      t(Y) %*% eta - sum((phi + YI) * log(phi +
                                            mu.i))
    return(logver)
  }

}

#' Maximum likelihood estimation
#'
#' @description Estimate parameters by quasi-Newton algorithms.
#' @param star Initial values for the parameters to be optimized over.
#' @param formula The structure matrix of covariates of dimension n x p (in models that include an intercept x
#' should contain a column of ones).
#' @param dataSet data
#' @param tab Logical. Print a summary of the coefficients,
#' standard errors and p-value for class "MNB".
#' @details Method "BFGS" is a quasi-Newton method, specifically that published
#' simultaneously in 1970 by Broyden, Fletcher, Goldfarb and Shanno. This uses function
#' values and gradients to build up a picture of the surface to be optimized.
#' @return Returns a list of summary statistics
#' of the fitted multivariate negative binomial model.
#' @author Jalmar M F Carrasco <carrascojalmar@gmail.com>,
#' Cristian M Villegas Lobos <master.villegas@gmail.com> and Lizandra C Fabio <lizandrafabio@gmail.com>
#' @references
#' \itemize{
#' \item Fabio, L., Paula, G. A., and de Castro, M. (2012). A Poisson mixed model
#' with nonormal random effect distribution. Computational Statistics and
#' Data Analysis, 56, 1499-1510.
#' \item Fabio, L. C., Villegas, C., Carrasco, J. M. F., and de Castro, M. (2023).
#' Diagnostic tools for a multivariate negative binomial
#' model for fitting correlated data with overdispersion. Communications in
#' Statistics - Theory and Methods, 52, 1833–1853.
#' \item Fabio, L. C., Villegas, C., Mamun, A. S., and Carrasco, J. M. F. (2025).
#' Residual analysis for discrete correlated data in the multivariate approach.
#' Brazilian Journal of Biometrics, 43, e43728.
#' }
#'
#' @examples
#'
#' \donttest{
#'
#' data(seizures)
#' head(seizures)
#'
#' star <-list(phi=1, beta0=1, beta1=1, beta2=1, beta3=1)
#'
#' mod1 <- fit.MNB(formula=Y ~ trt + period +
#' trt:period + offset(log(weeks)), star=star, dataSet=seizures)
#'
#' mod1
#'
#' seizures49 <- seizures[-c(241,242,243,244,245),]
#'
#' mod2 <- fit.MNB(formula=Y ~ trt + period +
#' trt:period + offset(log(weeks)), star=star, dataSet=seizures49)
#'
#' mod2
#'
#' }
#'
#' @export
#' @import stats

fit.MNB <- function(star, formula, dataSet, tab=TRUE) {
  op <- suppressWarnings(stats::optim(fn = l.MNB, gr = NULL, par = star,
                     method = "BFGS", hessian = TRUE,
                     formula = formula, dataSet = dataSet,
                     control = list(maxit = 500, trace = FALSE,
                                    fnscale = -1)))
  if(tab==FALSE){
    return(op)}
  else{
    se <- sqrt(diag(solve(-op$hessian)))
    z <- op$par/se
    pvalue <- 2 * (1 - stats::pnorm(abs(z)))
    TAB <- cbind(Estimate = op$par, Std.Error = se,
                 z.value = z, `Pr(>|z|)` = pvalue)
    mTab <- list( Lik = op$value,
                    Converged = op$convergence, Coefficients = TAB)
    return(mTab)
  }
}
