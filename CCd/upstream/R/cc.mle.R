cc.mle <- function(y) {

  n <- length(y)
  fun <- function(par, y, n) {
    mu <- par[1]  ;  lambda <- exp(par[2])
    - n * log( tanh(lambda * pi) ) - n * log(lambda) + sum( log( lambda^2 + (y - mu)^2) )
  }
  s <- 0.5 * abs( Rfast::nth(y, 3 * n/4) - Rfast::nth(y, n/4) )
  logs <- log(s)
  mod <- optim( c(mean(y), logs), y = y, n = n, fun, control = list(maxit = 5000) )
  param <- c( mod$par[1], exp(mod$par[2]) )
  names(param) <- c("mu", "lambda")
  list(param = param, loglik = -mod$value - n * log(pi))
}





