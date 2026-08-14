
### Windows CRAN needs too much time
if (.Platform$OS.type != "unix") quit()

library("mvtnorm")

set.seed(29)

chk <- function(...) stopifnot(isTRUE(all.equal(..., check.attributes = FALSE, tol = 1e-4)))

### N samples with N different covariance matrices

N <- 10
J <- 4

prm <- runif(N * J * (J + 1) / 2) + 1
m <- matrix(rnorm(N * J), nrow = J)
Z <- matrix(rnorm(N * J), ncol = N)
w <- matrix(runif((J - 1) * 100), nrow = J - 1)


thischeck <- expression({
  lt <- ltMatrices(matrix(prm[1: (N * J * (J + c(-1, 1)[dg + 1L]) / 2)], ncol = N), 
                   diag = dg)
  lt <- ltMatrices(lt, diag = dg, byrow = br)
  d <- Mult(lt, m)
  Y <- solve(lt, Z) + m
  lower <- Y - 2
  upper <- Y + 2

  l3 <- lpmvnorm(lower = lower, upper = upper, 
                  mean = m, invchol = lt, logLik = FALSE, w = w)
  l4 <- lpmvnorm(lower = lower, upper = upper,
                  mean = m, chol = solve(lt), logLik = FALSE, w = w)

  chk(l3, l4)

  l3d <- lpmvnorm(lower = lower, upper = upper,
                  invcholmean = d, invchol = lt, logLik = FALSE, w = w)
  l4d <- lpmvnorm(lower = lower, upper = upper,
                  invcholmean = d, chol = solve(lt), logLik = FALSE, w = w)

  chk(l3, l3d)
  chk(l4, l4d)

  ### check scores
  if (require("numDeriv", quietly = TRUE)) {

    f <- function(L) {
      L <- ltMatrices(L, diag = dg, byrow = br)
      lpmvnorm(lower = lower, upper = upper, mean = m, invchol = L, w = w)
    }

    s0 <- grad(f, unclass(lt))
    s1 <- slpmvnorm(lower = lower, upper = upper, mean = m, invchol = lt, w = w)

    chk(Lower_tri(ltMatrices(matrix(s0, ncol = N), diag = dg, byrow = br), diag = dg), 
        Lower_tri(s1$invchol, diag = dg))

    f <- function(L) {
      L <- ltMatrices(L, diag = dg, byrow = br)
      lpmvnorm(lower = lower, upper = upper, invcholmean = d, invchol = L, w = w)
    }

    s0 <- grad(f, unclass(lt))
    s1 <- slpmvnorm(lower = lower, upper = upper, invcholmean = d, invchol = lt, w = w)

    chk(Lower_tri(ltMatrices(matrix(s0, ncol = N), diag = dg, byrow = br), diag = dg), 
        Lower_tri(s1$invchol, diag = dg))

    f <- function(L) {
      L <- ltMatrices(L, diag = dg, byrow = br)
      lpmvnorm(lower = lower, upper = upper, mean = m, chol = L, w = w)
    }

    s0 <- grad(f, unclass(lt))
    s1 <- slpmvnorm(lower = lower, upper = upper, mean = m, chol = lt, w = w)

    chk(Lower_tri(ltMatrices(matrix(s0, ncol = N), diag = dg, byrow = br), diag = dg), 
        Lower_tri(s1$chol, diag = dg))



    f <- function(L) {
      L <- ltMatrices(L, diag = dg, byrow = br)
      lpmvnorm(lower = lower, upper = upper, invcholmean = d, chol = L, w = w)
    }

    s0 <- grad(f, unclass(lt))
    s1 <- slpmvnorm(lower = lower, upper = upper, invcholmean = d, chol = lt, w = w)

    chk(Lower_tri(ltMatrices(matrix(s0, ncol = N), diag = dg, byrow = br), diag = dg), 
        Lower_tri(s1$chol, diag = dg))

    f <- function(lwr)
      lpmvnorm(lower = lwr, upper = upper, mean = m, invchol = lt, w = w)

    s0 <- grad(f, lower)
    s1 <- slpmvnorm(lower = lower, upper = upper, mean = m, invchol = lt, w = w)

    chk(matrix(s0, ncol = N), s1$lower)

    f <- function(upr)
      lpmvnorm(lower = lower, upper = upr, mean = m, invchol = lt, w = w)

    s0 <- grad(f, upper)
    s1 <- slpmvnorm(lower = lower, upper = upper, mean = m, invchol = lt, w = w)

    chk(matrix(s0, ncol = N), s1$upper)


    f <- function(m)
      lpmvnorm(lower = lower, upper = upper, mean = m, invchol = lt, w = w)

    s0 <- grad(f, m)
    s1 <- slpmvnorm(lower = lower, upper = upper, mean = m, invchol = lt, w = w)

    chk(matrix(s0, ncol = N), s1$mean)

    f <- function(d)
      lpmvnorm(lower = lower, upper = upper, invcholmean = d, invchol = lt, w = w)

    s0 <- grad(f, d)
    s1 <- slpmvnorm(lower = lower, upper = upper, invcholmean = d, invchol = lt, w = w)

    chk(matrix(s0, ncol = N), s1$invcholmean)
  }
})

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
