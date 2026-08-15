#' @export
#' @import stats
pzibellt<-function (q, lambda, theta, pi, lower.tail = TRUE, log.p = FALSE) {
  if (any(lambda <= 0))
    stop(paste("lambda must be greater than 0"))
  if (any(theta <= 0))
    stop(paste("theta must be greater than 0"))
  if (any(pi < 0) | any(pi > 1))
    stop(paste("pi must be between 0 and 1"))
  if (any(q < 0))
    stop(paste("q must be 0 or greater than 0"))
  rval <- max(length(q), length(lambda), length(theta), length(pi))
  q <- rep(q, length = rval)
  pi <- rep(pi, length = rval)
  lambda <- rep(lambda, length = rval)
  theta <- rep(theta, length = rval)
  cdf <- rep(0, rval)
  cdf <- pbellt(q, lambda=lambda, theta=theta, lower.tail = TRUE, log.p = FALSE)
  cdf <- pi + (1 - pi) * cdf
  if (lower.tail == TRUE)
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE)
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}

