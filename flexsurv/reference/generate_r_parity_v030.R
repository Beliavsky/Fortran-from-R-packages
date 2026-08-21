## Reference calculations for flexsurv-fortran v0.3.0.
## Requires flexsurv; splines2 is used when available.

if (requireNamespace("splines2", quietly=TRUE)) {
  k <- c(0, 0.4, 1)
  x <- c(-0.5, 0.1, 0.37, 0.8, 1.5)
  b <- cbind(1, splines2::naturalSpline(
    x, knots=k[2], Boundary.knots=k[c(1,3)], intercept=FALSE))
  db <- cbind(0, splines2::naturalSpline(
    x, knots=k[2], Boundary.knots=k[c(1,3)], intercept=FALSE, derivs=1))
  print(list(splines2ns_basis=b, splines2ns_derivative=db))
}

if (requireNamespace("flexsurv", quietly=TRUE)) {
  tt <- c(0.2,0.4,0.7,1.0,1.3,1.7,2.1,2.9)
  fit <- flexsurv::flexsurvreg(survival::Surv(tt, rep(1,length(tt))) ~ 1, dist="exp")
  print(fit$res)
}

## Analytic shared-regression covariance reference used by test_shared_multistate.
theta <- c(log(0.1), 0.3)
V <- matrix(c(0.01,0.005,0.005,0.04),2,2)
t <- 2
H1 <- t*exp(theta[1]); H2 <- t*exp(sum(theta))
g1 <- c(H1,0); g2 <- c(H2,H2)
print(cross_cov=drop(t(g1)%*%V%*%g2))
