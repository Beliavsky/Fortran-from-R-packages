dcc <- function(y, mu = 0, lambda, logged = FALSE) {
  den <- log( tanh(lambda * pi) ) - log(pi) + log(lambda) - log( lambda^2 + (y - mu)^2 )
  if ( !logged )  den <- exp(den)
  den
}

