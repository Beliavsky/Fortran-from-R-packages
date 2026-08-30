
### Windows CRAN needs too much time
if (.Platform$OS.type != "unix") quit()

library("mvtnorm")

set.seed(29)

chk <- function(...) stopifnot(isTRUE(all.equal(..., check.attributes = FALSE, tol = 1e-4)))

### N samples with N different covariance matrices

N <- 2
J <- 8

prm <- runif(N * J * (J + 1) / 2) + 1
m <- matrix(rnorm(N * J), nrow = J)
Z <- matrix(rnorm(N * J), ncol = N)
W <- matrix(runif((J - 1) * 100), nrow = J - 1)

thischeck <- expression({
  lt <- ltMatrices(matrix(prm[1: (N * J * (J + c(-1, 1)[dg + 1L]) / 2)], ncol = N), 
                   diag = dg)
  lt <- ltMatrices(lt, diag = dg, byrow = br)
  d <- Mult(lt, m)
  Y <- solve(lt, Z) + m
  obs <- Y[idx1,]
  lower <- Y[idx2,] - 2
  upper <- Y[idx2,] + 2
  w <- W[seq_len(length(idx2) - 1),,drop = FALSE]

  objmL <- mvnorm(mean = m, invchol = lt)
  objmC <- mvnorm(mean = m, chol = solve(lt))

  l3 <- logLik(objmL, obs = obs, lower = lower, upper = upper, 
               logLik = FALSE, w = w)
  l4 <- logLik(objmC, obs = obs, lower = lower, upper = upper,
               logLik = FALSE, w = w)

  chk(l3, l4)

  objiL <- mvnorm(invcholmean = d, invchol = lt)
  objiC <- mvnorm(invcholmean = d, chol = solve(lt))

  l3d <- logLik(objiL, obs = obs, lower = lower, upper = upper,
                logLik = FALSE, w = w)
  l4d <- logLik(objiC, obs = obs, lower = lower, upper = upper,
                logLik = FALSE, w = w)

  chk(l3, l3d)
  chk(l4, l4d)

  ### check scores
  if (require("numDeriv", quietly = TRUE)) {

    f <- function(L) {
      L <- ltMatrices(L, diag = dg, byrow = br)
      obj <- mvnorm(mean = m, invchol = L)
      logLik(obj, obs = obs, lower = lower, upper = upper, w = w)
    }

    s0 <- grad(f, unclass(lt))
    s1 <- lLgrad(objmL, obs = obs, lower = lower, upper = upper, w = w)

    chk(Lower_tri(ltMatrices(matrix(s0, ncol = N), diag = dg, byrow = br), diag = dg), 
        Lower_tri(s1$scale, diag = dg))

    f <- function(L) {
      L <- ltMatrices(L, diag = dg, byrow = br)
      obj <- mvnorm(invcholmean = d, invchol = L)
      logLik(obj, obs = obs, lower = lower, upper = upper, w = w)
    }

    s0 <- grad(f, unclass(lt))
    s1 <- lLgrad(objiL, obs = obs, lower = lower, upper = upper, w = w)

    chk(Lower_tri(ltMatrices(matrix(s0, ncol = N), diag = dg, byrow = br), diag = dg), 
        Lower_tri(s1$scale, diag = dg))

    f <- function(L) {
      L <- ltMatrices(L, diag = dg, byrow = br)
      obj <- mvnorm(mean = m, chol = L)
      logLik(obj, obs = obs, lower = lower, upper = upper, w = w)
    }

    s0 <- grad(f, unclass(solve(lt)))
    s1 <- lLgrad(objmC, obs = obs, lower = lower, upper = upper, w = w)

    chk(Lower_tri(ltMatrices(matrix(s0, ncol = N), diag = dg, byrow = br), diag = dg), 
        Lower_tri(s1$scale, diag = dg))

    f <- function(L) {
      L <- ltMatrices(L, diag = dg, byrow = br)
      obj <- mvnorm(invcholmean = d, chol = L)
      logLik(obj, obs = obs, lower = lower, upper = upper, w = w)
    }

    s0 <- grad(f, unclass(solve(lt)))
    s1 <- lLgrad(objiC, obs = obs, lower = lower, upper = upper, w = w)

    chk(Lower_tri(ltMatrices(matrix(s0, ncol = N), diag = dg, byrow = br), diag = dg), 
        Lower_tri(s1$scale, diag = dg))

    f <- function(x)
      logLik(objmL, obs = x, lower = lower, upper = upper, w = w)

    s0 <- grad(f, obs)
    s1 <- lLgrad(objmL, obs = obs, lower = lower, upper = upper, w = w)

    chk(matrix(s0, ncol = N), s1$obs)

    f <- function(x)
      logLik(objiL, obs = x, lower = lower, upper = upper, w = w)

    s0 <- grad(f, obs)
    s1 <- lLgrad(objiL, obs = obs, lower = lower, upper = upper, w = w)

    chk(matrix(s0, ncol = N), s1$obs)

    f <- function(lwr)
      logLik(objmL, obs = obs, lower = lwr, upper = upper, w = w)

    s0 <- grad(f, lower)
    s1 <- lLgrad(objmL, obs = obs, lower = lower, upper = upper, w = w)

    chk(matrix(s0, ncol = N), s1$lower)

    f <- function(lwr)
      logLik(objiL, obs = obs, lower = lwr, upper = upper, w = w)

    s0 <- grad(f, lower)
    s1 <- lLgrad(objiL, obs = obs, lower = lower, upper = upper, w = w)

    chk(matrix(s0, ncol = N), s1$lower)

    f <- function(upr)
      logLik(objmL, obs = obs, lower = lower, upper = upr, w = w)

    s0 <- grad(f, upper)
    s1 <- lLgrad(objmL, obs = obs, lower = lower, upper = upper, w = w)

    chk(matrix(s0, ncol = N), s1$upper)

    f <- function(upr)
      logLik(objiL, obs = obs, lower = lower, upper = upr, w = w)

    s0 <- grad(f, upper)
    s1 <- lLgrad(objiL, obs = obs, lower = lower, upper = upper, w = w)

    chk(matrix(s0, ncol = N), s1$upper)

    f <- function(m) {
      obj <- mvnorm(mean = m, invchol = lt)
      logLik(obj, obs = obs, lower = lower, upper = upper, w = w)
    }

    s0 <- grad(f, m)
    s1 <- lLgrad(objmL, obs = obs, lower = lower, upper = upper, w = w)

    chk(matrix(s0, ncol = N), s1$mean)

    f <- function(d) {
      obj <- mvnorm(invcholmean = d, invchol = lt)
      logLik(obj, obs = obs, lower = lower, upper = upper, w = w)
    }

    s0 <- grad(f, d)
    s1 <- lLgrad(objiL, obs = obs, lower = lower, upper = upper, w = w)

    chk(matrix(s0, ncol = N), s1$invcholmean)
  }
})


idx <- seq_len(J)
idx1 <- idx[1:4]
idx2 <- idx[-(1:4)]

dg <- TRUE
br <- FALSE
eval(thischeck)

dg <- FALSE
br <- FALSE
eval(thischeck)

dg <- FALSE
br <- TRUE
eval(thischeck)

dg <- FALSE
br <- FALSE
eval(thischeck)

idx1 <- idx[-(1:4)]
idx2 <- idx[1:4]

dg <- TRUE
br <- FALSE
eval(thischeck)

dg <- FALSE
br <- FALSE
eval(thischeck)

dg <- FALSE
br <- TRUE
eval(thischeck)

dg <- FALSE
br <- FALSE
eval(thischeck)

