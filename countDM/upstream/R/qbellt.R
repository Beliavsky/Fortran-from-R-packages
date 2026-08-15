#' @export
#' @import stats

qbellt<-function (p, lambda, theta,lower.tail = TRUE, log.p = FALSE) {
  if (any(lambda <= 0))
    stop(paste("lambda must be greater than 0"))
  if (any(theta <= 0))
    stop(paste("theta must be greater than 0"))
  if (any(p < 0) | any(p > 1))
    warning("pi must be in 0, 1")

  n <- length(p)
  k <- rep(0, n)
  for ( i in 1:n ) {
    q <- 0
    while ( q <= p[i] ) {
      q <- pbellt(k[i], lambda, theta)
      k[i] <- k[i] + 1
    }
  }
  if ( log.p ) {
    k <- log(k - 1)
  } else  k <- k - 1

  k
}
