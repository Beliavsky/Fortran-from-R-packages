
library("mvtnorm")

set.seed(29)

chk <- function(...) stopifnot(isTRUE(all.equal(..., check.attributes = FALSE)))

### N samples with N different covariance matrices

N <- 10
J <- 5
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
  ll1 <- sum(dnorm(Mult(lt, Y - m), log = TRUE)) + sum(log(diagonals(lt)))

  S <- as.array(Tcrossprod(solve(lt)))
  ll2 <- sum(l2 <- sapply(1:N, function(i) mvtnorm:::dmvnorm(x = Y[,i], mean = m[,i], sigma = S[,,i], log = TRUE)))
  chk(ll1, ll2)

  l3 <- ldmvnorm(obs = Y, mean = m, invchol = lt, logLik = FALSE)
  l4 <- ldmvnorm(obs = Y, mean = m, chol = solve(lt), logLik = FALSE)

  chk(l2, l3)
  chk(l2, l4)

  l3d <- ldmvnorm(obs = Y, invcholmean = d, invchol = lt, logLik = FALSE)
  l4d <- ldmvnorm(obs = Y, invcholmean = d, chol = solve(lt), logLik = FALSE)

  chk(l3, l3d)
  chk(l4, l4d)

  ll1 <- sum(dnorm(Mult(lt[1,], Y - m[,1]), log = TRUE)) + N * sum(log(diagonals(lt[1,])))

  S <- as.array(Tcrossprod(solve(lt)))
  ll2 <- sum(l2 <- sapply(1:N, function(i) mvtnorm:::dmvnorm(x = Y[,i], mean = m[,1], sigma = S[,,1], log = TRUE)))
  chk(ll1, ll2)

  l3 <- ldmvnorm(obs = Y, mean = m[,1], invchol = lt[1,], logLik = FALSE)
  l4 <- ldmvnorm(obs = Y, mean = m[,1], chol = solve(lt[1,]), logLik = FALSE)

  chk(l2, l3)
  chk(l2, l4)

  l3d <- ldmvnorm(obs = Y, invcholmean = d[,1], invchol = lt[1,], logLik = FALSE)
  l4d <- ldmvnorm(obs = Y, invcholmean = d[,1], chol = solve(lt[1,]), logLik = FALSE)

  chk(l3, l3d)
  chk(l4, l4d)

  ### check scores
  if (require("numDeriv", quietly = TRUE)) {

    f <- function(L) {
      L <- ltMatrices(L, diag = dg, byrow = br)
      ldmvnorm(obs = Y, mean = m, invchol = L)
    }

    s0 <- grad(f, unclass(lt))
    s1 <- sldmvnorm(obs = Y, mean = m, invchol = lt)

    chk(Lower_tri(ltMatrices(matrix(s0, ncol = N), diag = dg, byrow = br), diag = dg), 
        Lower_tri(s1$invchol, diag = dg))

    f <- function(L) {
      L <- ltMatrices(L, diag = dg, byrow = br)
      ldmvnorm(obs = Y, invcholmean = d, invchol = L)
    }

    s0 <- grad(f, unclass(lt))
    s1 <- sldmvnorm(obs = Y, invcholmean = d, invchol = lt)

    chk(Lower_tri(ltMatrices(matrix(s0, ncol = N), diag = dg, byrow = br), diag = dg), 
        Lower_tri(s1$invchol, diag = dg))

    f <- function(L) {
      L <- ltMatrices(L, diag = dg, byrow = br)
      ldmvnorm(obs = Y, mean = m, chol = L)
    }

    s0 <- grad(f, unclass(lt))
    s1 <- sldmvnorm(obs = Y, mean = m, chol = lt)

    chk(Lower_tri(ltMatrices(matrix(s0, ncol = N), diag = dg, byrow = br), diag = dg), 
        Lower_tri(s1$chol, diag = dg))


    f <- function(L) {
      L <- ltMatrices(L, diag = dg, byrow = br)
      ldmvnorm(obs = Y, invcholmean = d, chol = L)
    }

    s0 <- grad(f, unclass(lt))
    s1 <- sldmvnorm(obs = Y, invcholmean = d, chol = lt)

    chk(Lower_tri(ltMatrices(matrix(s0, ncol = N), diag = dg, byrow = br), diag = dg), 
        Lower_tri(s1$chol, diag = dg))

    f <- function(x)
      ldmvnorm(obs = x, mean = m, invchol = lt)

    s0 <- grad(f, Y)
    s1 <- sldmvnorm(obs = Y, mean = m, invchol = lt)

    chk(matrix(s0, ncol = N), s1$obs)

    f <- function(x)
      ldmvnorm(obs = x, invcholmean = d, invchol = lt)

    s0 <- grad(f, Y)
    s1 <- sldmvnorm(obs = Y, invcholmean = d, invchol = lt)

    chk(matrix(s0, ncol = N), s1$obs)

    f <- function(m)
      ldmvnorm(obs = Y, mean = m, invchol = lt)

    s0 <- grad(f, m)
    s1 <- sldmvnorm(obs = Y, mean = m, invchol = lt)

    chk(matrix(s0, ncol = N), -s1$obs)

    f <- function(d)
      ldmvnorm(obs = Y, invcholmean = d, invchol = lt)

    s0 <- grad(f, d)
    s1 <- sldmvnorm(obs = Y, invcholmean = d, invchol = lt)

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

