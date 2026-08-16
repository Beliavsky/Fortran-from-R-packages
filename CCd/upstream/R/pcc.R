pcc <- function(y, mu = 0, lambda) {

  fun <- function(y, mu, lambda) {
    if ( y == 0 ) {
      res <- 0.5 + 0.5 * CCd::dcc(0, mu, lambda)
    } else if ( y < 0 ) {
      res <- 0.5 - 0.5 * sum( CCd::dcc( (y + 1):(abs(y) - 1), mu, lambda ) )
    } else if (y > 0 ) {
      res <- 0.5 + 0.5 * CCd::dcc(0, mu, lambda) + sum( CCd::dcc(1:y, mu, lambda) )
    }
    res
  }

  n <- length(y)  ;   p <- numeric(n)
  for (i in 1:n)  p[i] <- fun(y[i], mu, lambda)
  p
}

