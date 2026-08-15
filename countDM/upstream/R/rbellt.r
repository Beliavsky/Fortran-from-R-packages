#' @export
#' @import stats
rbellt<-function (n, lambda, theta) {
  if (any(lambda <= 0))
    stop(paste("lambda must greated than 0"))
  if (any(theta <= 0))
    stop(paste("theta must greated than 0"))
  if (any(n <= 0))
    stop(paste("n must be a positive integer"))
  n <- ceiling(n)
  p <- runif(n)
  z <- qbellt(p, lambda=lambda, theta=theta, lower.tail = TRUE, log.p = FALSE)
  as.integer(z)
}
