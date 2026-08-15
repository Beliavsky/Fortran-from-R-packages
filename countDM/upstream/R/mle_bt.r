#' @export
#' @import maxLik
#' @import stats
#' @import numbers
#' @importFrom miscTools stdEr
mle_bt <- function(x,lambda,theta) {
  n <- length(x)

loglik<-function(p){
    lambda <- p[1]
    theta <- p[2]

loglik <- sum(x) * log(p[1]) + n * p[2] * (1 - exp(p[1])) - sum(lgamma(x + 1)) + sum(log( sapply(x, theta=p[2], TP)))
}

res <- suppressWarnings(maxLik::maxLik(loglik,start=c(lambda,theta)))
  aux <- cbind(coef(res), miscTools::stdEr(res))
  colnames(aux) = c("MLE", "SE")

  aux1 = cbind(AIC(res))
  colnames(aux1) = c("AIC")
  list(Estimates = aux, AIC= aux1)
}
