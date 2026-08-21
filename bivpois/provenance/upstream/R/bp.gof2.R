########################
### Goodness of fit for the bivariate poisson
### Vectorised version
####  3/2015
#### mtsagris@yahoo.gr
########################

bp.gof2 <- function(x1, x2 = NULL, R = 999) {
  if ( is.null(x2) ) {
    x2 <- x1[, 2]
    x1 <- x1[, 1]
  }
  x1 <- as.numeric(x1)   ;    x2 <- as.numeric(x2)
  ## x1 and x2 are the two variables
  runtime <- proc.time()
  n <- length(x1)  ## sample size
  r <- cor(x1, x2)  ## Pearson correlation coefficient
  m1 <- mean(x1)     ;    m2 <- mean(x2)
  Ib <- n/(1 - r^2) * ( var(x1) / m1 + var(x2) / m2 - 2 * r^2 * sqrt( var(x1) / m1 * var(x2) / m2 ) )  ## test statistic
  tab <- table(x1, x2)

  lambda <- bivpois::bp.mle(x1, x2)$lambda
  z3 <- matrix( rpois(R * n, lambda[3]), ncol = R )
  z1 <- matrix( rpois(R * n, lambda[1]), ncol = R ) + z3
  z2 <- matrix( rpois(R * n, lambda[2]), ncol = R ) + z3
  m1 <- Rfast::colmeans(z1)     ;    m2 <- Rfast::colmeans(z2)
  v1 <- Rfast::colVars(z1)      ;    v2 <- Rfast::colVars(z2)
  sxy <- Rfast::colsums(z1 * z2)
  rb <- (sxy - n * m1 * m2) / ( (n - 1) * sqrt( v1 * v2 ) )
  tb <- n/(1 - rb^2) * ( v1 / m1 + v2 / m2 - 2 * rb^2 * sqrt( v1 / m1 * v2 / m2 ) )

  pvalue <- ( sum(tb > Ib) + 1 ) / (R + 1)
  runtime <- proc.time() - runtime
  list(runtime = runtime, tab = tab, pvalue = pvalue)
}
