#' @export
#' @import stats
rzibellt<-function (n, lambda, theta, pi) {
  if (any(lambda <= 0))
    stop(paste("lambda must greated than 0"))
  if (any(theta <= 0))
    stop(paste("theta must greated than 0"))
  if (any(pi < 0) | any(pi > 1))
    warning("pi must be in 0, 1")
  if (any(n <= 0))
    stop(paste("n must be a positive integer"))
  n <- ceiling(n)
  p <- runif(n)
  z <- qzibellt(p, lambda=lambda, theta=theta, pi=pi, lower.tail = TRUE, log.p = FALSE)
  as.integer(z)
}
