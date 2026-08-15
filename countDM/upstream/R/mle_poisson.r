#' @export
#' @import maxLik
#' @import stats
#' @importFrom miscTools stdEr

mle_poisson <- function(x, theta) {


n <- length(x)

  loglik<-function(p){

    theta <- p[1]
    loglik <- sum(x)*log(p[1])-n*p[1]-sum(lgamma(x+1))

  }
  res <- suppressWarnings(maxLik(loglik,start=c(theta)))
  aux <- cbind(coef(res), miscTools::stdEr(res))
  colnames(aux) = c("MLE", "SE")

  aux1 = cbind(AIC(res))
  colnames(aux1) = c("AIC")
  list(Estimates = aux, AIC= aux1)
}
