bp.contour <- function(x1, x2 = NULL, lambda) {
  if ( is.null(x2) ) {
    x2 <- x1[, 2]
    x1 <- x1[, 1]
  }
  x1 <- as.numeric(x1)   ;    x2 <- as.numeric(x2)
  ## x1 and x2 are the two variables
  ## lambda contains the three values of the three parameters
  lam1 <- lambda[1]  ;  lam2 <- lambda[2]   ;   lam3 <- lambda[3]
  z1 <- seq(max(min(x1) - 3, 0), max(x1) + 3)
  n1 <- length(z1)
  z2 <- seq(max(min(x2) - 3, 0), max(x2) + 3)
  n2 <- length(z2)
  mat <- matrix(nrow = n1, ncol = n2)
  l1 <- log(lam1)  ;  l2 <- log(lam2)
  ls <-  -(lam1 + lam2 + lam3)
  rho <- lam3 / sqrt(lam1 * lam2)
  for (i in 1:n1) {
    for (j in 1:n2) {
      f1 <- ls + z1[i] * l1 - lgamma(z1[i] + 1) + z2[j] * l2 - lgamma(z2[j] + 1)
      k <- 0:min(z1[i], z2[j])
      f2 <- log( sum( choose(z1[i], k) * choose(z2[j], k) * factorial(k) * rho^k ) )
      f <- f1 + f2
      mat[i, j] <- exp(f)
    }
  }
  contour(z1, z2, mat, nlevels = 10, col = 2, xlab = "X", ylab = "Y")
  points(x1, x2, pch = 20)
}
