########################
### Goodness of fit for the bivariate poisson
####  3/2015
#### mtsagris@yahoo.gr
########################

bp.gof <- function(x1, x2 = NULL, R = 999) {
  if ( is.null(x2) ) {
    x2 <- x1[, 2]
    x1 <- x1[, 1]
  }
  x1 <- as.numeric(x1)   ;    x2 <- as.numeric(x2)
  ## x1 and x2 are the two variables
  runtime <- proc.time()
  n <- length(x1)  ## sample size
  r <- cor(x1, x2)  ## Pearson correlation coefficient
  m1 <- sum(x1) / n     ;     m2 <- sum(x2)/n
  v1 <- ( sum(x1^2) - n * m1^2 ) / (n - 1)
  v2 <- ( sum(x2^2) - n * m2^2 ) / (n - 1)

  Ib <- n/(1 - r^2) * ( v1 / m1 + v2 / m2 - 2 * r^2 * sqrt(v1 / m1 * v2 / m2) )  ## test statistic
  tab <- table(x1, x2)
  tb <- numeric(R)
  lambda <- bivpois::bp.mle(x1, x2)$lambda

  for (i in 1:R) {
    z3 <- rpois(n, lambda[3])
    z1 <- rpois(n, lambda[1]) + z3
    z2 <- rpois(n, lambda[2]) + z3
    r <- cor(z1, z2)
    m1 <- sum(z1)/n    ;    m2 <- sum(z2)/n
    s1 <- ( sum(z1^2) - n * m1^2 ) / (n - 1)
    s2 <- ( sum(z2^2) - n * m2^2 ) / (n - 1)
    tb[i] <- n/(1 - r^2) * ( s1 / m1 + s2 / m2 - 2 * r^2 * sqrt( s1 / m1 * s2 / m2 ) )
  }

  pvalue <- (sum(tb > Ib) + 1)/(R + 1)
  runtime <- proc.time() - runtime
  list(runtime = runtime, tab = tab, pvalue = pvalue)
}
