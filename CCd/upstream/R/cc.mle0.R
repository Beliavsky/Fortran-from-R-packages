cc.mle0 <- function(y, tol = 1e-7) {

  n <- length(y)
  fun <- function(lambda, y, n) {
    n * log( tanh(lambda * pi) ) + n * log(lambda) - sum( log( lambda^2 + y^2) )
  }
  mod <- optimize(fun, c(0, 1000), y = y, n = n, maximum = TRUE, tol = tol )
  list(lambda = mod$maximum, loglik = mod$objective - n * log(pi))
}


