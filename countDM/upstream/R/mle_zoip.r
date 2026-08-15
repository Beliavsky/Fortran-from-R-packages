#' @export
#' @import maxLik
#' @import stats
#' @importFrom miscTools stdEr

mle_zoip <- function(x, alpha, beta, theta) {

no <- sum(x == 0)
noo <- sum(x == 1)
n <- length(x)
n1 <- (n - no - noo)
x1 <- x[x > 1]
  loglik<-function(p){

    alpha <- p[1]
	beta <- p[2]
	theta <- p[3]
    loglik <- no * log(p[1] + (1 - p[1]-p[2]) * exp(-p[3])) + noo * log(p[2] + (1 - p[1]- p[2]) *p[3]*exp(-p[3]))+ n1 * log(1-p[1]- p[2]) + sum(dpois(x1, p[3], log = TRUE))

  }
  res <- suppressWarnings(maxLik(loglik,start=c(alpha, beta, theta)))
  aux <- cbind(coef(res), miscTools::stdEr(res))
  colnames(aux) = c("MLE", "SE")

  aux1 = cbind(AIC(res))
  colnames(aux1) = c("AIC")
  list(Estimates = aux, AIC= aux1)
}
