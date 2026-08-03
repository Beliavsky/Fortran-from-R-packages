#### From Yanhao Shi's thesis -- R/Functions/CLA_analysis.R

## Vectorized in first argument by Martin Maechler

### Version 1 -- of the code, simply "vectorize via lapply" -- not so fast

findMu <- function(Sig0, result, covar, tol.unir = 1e-6, equal.tol = 1e-6) {
    R <- lapply(Sig0, findMu.1,
                result=result, covar=covar, tol.unir=tol.unir, equal.tol=equal.tol)
    ## each R[[i]] has (Mu, weight)
    list(Mu     = vapply(R, `[[`, numeric(1),           "Mu"),
         weight = vapply(R, `[[`, numeric(nrow(covar)), "weight"))
}

findSig <- function(Mu0, result, covar, equal.tol = 1e-6) { # equal.tol > 1e-16
    R <- lapply(Mu0, findSig.1,
                result=result, covar=covar, equal.tol=equal.tol)
    ## each R[[i]] has (Sig, weight):
    list(Sig    = vapply(R, `[[`, numeric(1),           "Sig"),
         weight = vapply(R, `[[`, numeric(nrow(covar)), "weight"))
}


##' Now simplified and very close to findSig0() -- just returning weights additionally
findSig.1 <- function(Mu0, result, covar, equal.tol) { # equal.tol > 1e-16
  stopifnot(length(Mu0) == 1L)
  ms.w <- result$MS_weight
  nd <- order(ms.w[, "Mu"])
  n <- length(nd)
  ## revert order for all (w, mu, sig):
  ms.w <- ms.w[nd, ]
  w <- result$weights_set[, nd]
  mu.w  <- ms.w[, "Mu"]
  sig.w <- ms.w[, "Sig"]
  if(Mu0 < min(mu.w) || Mu0 > max(mu.w))
      stop(sprintf("Mu0 must be in [%g, %g]", min(mu.w), max(mu.w)))
  i <- findInterval(Mu0, mu.w) # [..)...[..]
  if(i == n || # Mu0 is max(mu.w)
     isTRUE(all.equal(mu.w[i], mu.w[i+1], tol = equal.tol))) {
    list(Sig = sig.w[i], weight = w[, i])
  } else {
    a <- (Mu0 - mu.w[i+1])/(mu.w[i] - mu.w[i+1])
    # solve for a : Mu0 = a* Mu1 + (1-a)* Mu2 :
    wm <- a*w[, i] + (1-a)*w[, i+1]
    list(Sig = sqrt(colSums(wm *(covar %*% wm) )),
         weight = wm)
  }
}

##' This is much cleaner and logical than the original findSig() ==> use a version of it
##' A version of findSig() only used inside  findMu():
findSig0 <- function(Mu0, result, covar, equal.tol = 1e-8){ # equal.tol > 1e-16
  stopifnot(length(Mu0) == 1L)
  ms.w <- result$MS_weight
  nd <- order(ms.w[, "Mu"])
  n <- length(nd)
  ## revert order for all (w, mu, sig):
  ms.w <- ms.w[nd, ]
  w <- result$weights_set[, nd]
  mu.w  <- ms.w[, "Mu"]
  sig.w <- ms.w[, "Sig"]
  i <- findInterval(Mu0, mu.w) # [..)...[..]
  if(i == n || isTRUE(all.equal(mu.w[i], mu.w[i+1], tol = equal.tol))) {
    sig.w[i]
  } else {
    a <- (Mu0 - mu.w[i+1])/(mu.w[i] - mu.w[i+1])
    # solve for a : Mu0 = a* Mu1 + (1-a)* Mu2 :
    w0 <- a*w[, i] + (1-a)*w[, i+1]
    sqrt(colSums(w0 *(covar %*% w0) )) #return Sig0
  }
}

findMu.1 <- function(Sig0, result, covar, tol.unir = 1e-6, equal.tol = 1e-6) {
  stopifnot(length(Sig0) == 1L)
  ms.w <- result$MS_weight
  nd <- order(ms.w[, "Sig"])
  n <- length(nd)
  mu.w  <- ms.w[nd, "Mu"]
  sig.w <- ms.w[nd, "Sig"]
  w <- result$weights_set[, nd]
  if(Sig0 < min(sig.w) || Sig0 > max(sig.w))
      stop(sprintf("Sig0 must be in [%g, %g]", min(sig.w), max(sig.w)))
  i <- findInterval(Sig0, sig.w)
  if(i == n || # Mu0 is max(mu.w)
     isTRUE(all.equal(mu.w[i], mu.w[i+1], tol = equal.tol))) {
    list(Mu = mu.w[i], weight = w[, i])
  } else { ## FIXME: here are using default equal.tol = 1e-8 in findSig0() !
    r <- uniroot(function(mu) findSig0(mu, result, covar) - Sig0,
                 interval = mu.w[c(i, i+1)], tol=tol.unir)
    Mu0 <- r$root
    # solve for a : mu = a* Mu1 + (1-a)* Mu2 :
    a <- (Mu0 - mu.w[i+1])/(mu.w[i] - mu.w[i+1])
    w0 <- a* w[, i] + (1-a)*w[, i+1]
    list(Mu = Mu0, weight = w0)
  }
}


### Version 2 --- vectorize inside -- only those parts that need it

## not yet !!!  === For first part, look at
## ===========      findSigMu.R+
##                  ============
if(FALSE) ## will be .... not yet
findMu <- function(Sig0, result, covar, tol.unir = 1e-6, equal.tol = 1e-6){
  ms.w <- result$MS_weight
  nd <- order(ms.w[, "Sig"])
  n <- length(nd)
   mu.w <- ms.w[nd, "Mu"]
  sig.w <- ms.w[nd, "Sig"]
  w <- result$weights_set[, nd]
  if(Sig0 < min(sig.w) || Sig0 > max(sig.w))
      stop(sprintf("Sig0 must be in [%g, %g]", min(sig.w), max(sig.w)))
  ini <- findInterval(Sig0, sig.w)
  m <- length(Sig0)
  mu0 <- numeric(m)
  wt0 <- matrix(NA_real_, n, m)
  iBnd <- vapply(ini, function(i) {
      (i == n || # Mu0 is max(mu.w)
       isTRUE(all.equal(mu.w[i], mu.w[i+1], tol = equal.tol))) # duplicate turning pt
  }, NA)
  if(any(iBnd)) {
      i <- ini[iBnd]
      mu0[ iBnd] <- mu.w[i]
      wt0[,iBnd] <- w[, i]
  }
  if(any(iIn <- !iBnd)) { # regular case
      i <- ini[iIn]
      mus <- vapply(..., function(i) {
          r <- uniroot(function(mu) findSig0(mu, result, covar) - Sig0[i],
                       interval = mu.w[c(i, i+1)], tol = tol.unir)
          ## TODO check convergence?
          r$root
      }, numeric(1))

      ## solve for a : mu = a* Mu1 + (1-a)* Mu2  ==>
      a <- (Mu0 - mu.w[i+1])/(mu.w[i] - mu.w[i+1])
      w0 <- a* w[, i] + (1-a)*w[, i+1]
      mu0[ii] <- mu.w[i]
      wt0[,ii] <- w[, i]

  }
  list(Mu = Mu0, weight = w0)
}










