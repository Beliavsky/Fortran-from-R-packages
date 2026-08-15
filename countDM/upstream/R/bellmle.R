#' @export
#' @import lamW
#' @import stats
#' @import numbers
bell_mle <- function(x) {

  n <- length(x)

  sx <- sum(x)
  theta <- lamW::lambertW0(sx/n)


  loglik <- sx*log(theta)+n*(1-exp(theta))+sum(log(sapply(x,numbers::bell)))-sum(lgamma(x + 1))

  list(loglik =loglik, theta = theta)

}
