qcc <- function(p, mu, lambda) {

  fun <- function(y, p, mu, lambda) {
    abs( p - CCd::pcc(y, mu, lambda) )
  }

  n <- length(p)
  quan <- numeric(n)
  if (n == 1) {
    quan <- ceiling( optimize(fun, c(-1e+6, 1e+6), p = p, mu = mu, lambda = lambda, tol = 1e-12)$minimum )
  } else {
    for (i in 1:n)
      quan[i] <- ceiling( optimize(fun, c(-1e+6, 1e+6), p = p[i], mu = mu, lambda = lambda, tol = 1e-12)$minimum )
  }
  quan
}
