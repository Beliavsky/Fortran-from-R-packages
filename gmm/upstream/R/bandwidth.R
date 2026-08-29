#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


bwWilhelm <- function(x, order.by = NULL, kernel = c("Quadratic Spectral",
  "Bartlett", "Parzen", "Tukey-Hanning"), approx = c("AR(1)", "ARMA(1,1)"),
  weights = NULL, prewhite = 1, ar.method = "ols", data=list()) 
{
  ## ensure that NAs are omitted
  if(is.list(x) && !is.null(x$na.action)) class(x$na.action) <- "omit"

  # ensure that x is a gmm object
  stopifnot(attr(x,"class") == "gmm")

  kernel <- match.arg(kernel)
  approx <- match.arg(approx)
  prewhite <- as.integer(prewhite)

  # the following line applies only to the case class(res)=="gmm"
  umat <- x$gt
  n <- nrow(umat)
  k <- ncol(umat)

  # the following line applies only to the case class(res)=="gmm"
  G <- x$G
  p <- ncol(G)

  # return Andrews' bandwidth in the just-identified case
  
  if (p==k)
      {
          class(umat) <- "gmmFct"
          return(bwAndrews(umat, order.by=order.by, kernel=kernel,
                           approx=approx, weights=weights, prewhite=prewhite,
                           ar.method=ar.method, data=data))
      }

  if(!is.null(order.by))
  {
    if(inherits(order.by, "formula")) {
      z <- model.matrix(order.by, data = data)
      z <- as.vector(z[,ncol(z)])
    } else {
      z <- order.by
    }
    index <- order(z)
  } else {
    index <- 1:n
  }

  umat <- umat[index, , drop = FALSE]

  ## compute weights (w_a, see p. 834)
  ## (try to set the intercept weight to 0)
  #### could be ignored by using: weights = 1
  
  if(is.null(weights)) {
    weights <- rep(1, p)
    unames <- colnames(umat)
    if(!is.null(unames) && "(Intercept)" %in% unames) 
      weights[which(unames == "(Intercept)")] <- 0
    else {
      res <- try(as.vector(rowMeans(estfun(x)/model.matrix(x), na.rm = TRUE)), silent = TRUE)
      if(inherits(res, "try-error")) res <- try(residuals(x), silent = TRUE)
      if(!inherits(res, "try-error")) weights[which(colSums((umat - res)^2) < 1e-16)] <- 0
    }
    if(isTRUE(all.equal(weights, rep(0, p)))) weights <- rep(1, p)
  } else {
    weights <- rep(weights, length.out = p)
  }
  if(length(weights) < 2) weights <- 1

  ## pre-whitening
  if(prewhite > 0) {
    var.fit <- ar(umat, order.max = prewhite, demean = FALSE, aic = FALSE, method = ar.method)
    if(inherits(var.fit, "try-error")) stop(sprintf("VAR(%i) prewhitening of estimating functions failed", prewhite))
    umat <- as.matrix(na.omit(var.fit$resid))
    n <- n - prewhite
  }
 
  # define kernel constants
  kernConst <- switch(kernel,
    "Bartlett"           = list(q = 1, g_q = 1,        mu1 = 1,       mu2 = 2/3),
    "Parzen"             = list(q = 2, g_q = 6,        mu1 = 0.75,    mu2 = 0.539286),
    "Tukey-Hanning"      = list(q = 2, g_q = pi^2/4,   mu1 = 1,       mu2 = 0.75),
    "Quadratic Spectral" = list(q = 2, g_q = 1.421223, mu1 = 1.25003, mu2 = 0.999985))

  # fit approximating model to moments
  if(approx == "AR(1)") {

    fitAR1 <- function(x) {
      rval <- ar(x, order.max = 1, aic = FALSE, method = "ols")
      rval <- c(rval$ar, sqrt(rval$var.pred))
      names(rval) <- c("rho", "sigma")
      return(rval)
    }
    
    ar.coef <- apply(umat, 2, fitAR1)

    Omega0 <- ar.coef["sigma", ]^2 / (1-ar.coef["rho", ])^2
    Omega_q <- if(kernConst$q == 1) {
      diag(2 * ar.coef["sigma", ]^2 * ar.coef["rho", ] / ((1 - ar.coef["rho", ])^3 * (1 + ar.coef["rho", ])))
    } else {
      diag(2 * ar.coef["sigma", ]^2 * ar.coef["rho", ] / (1 - ar.coef["rho", ])^4) 
    }

  } else {

    fitARMA11 <- function(x) {
      rval <- arima(x, order = c(1, 0, 1), include.mean = FALSE)
      rval <- c(rval$coef, sqrt(rval$sigma2))
      names(rval) <- c("rho", "psi", "sigma")
      return(rval)
    }
    
    arma.coef <- apply(umat, 2, fitARMA11)

    Omega0 <- (1 + ar.coef["psi", ])^2 * ar.coef["sigma", ]^2 / (1 - ar.coef["rho", ])^2
    Omega_q <- if(kernConst$q == 1) {
      diag(2 * (1 + ar.coef["psi", ] * ar.coef["rho", ]) * (ar.coef["psi", ] + ar.coef["rho", ]) * ar.coef["sigma", ]^2 / ((1 - ar.coef["rho", ])^3 * (1 + ar.coef["rho", ])))
    } else {
      diag(2 * (1 + ar.coef["psi", ] * ar.coef["rho", ]) * (ar.coef["psi", ] + ar.coef["rho", ]) * ar.coef["sigma", ]^2 / (1 - ar.coef["rho", ])^4)
    }

  }

  # compute remaining bandwidth components
  W <- diag(weights)
  Omega0inv <- diag(1/Omega0)
  Sigma0 <- solve( crossprod(Omega0inv%*%G, G) )
  H0 <- Sigma0 %*% t(G) %*% Omega0inv
  P0 <- Omega0inv - Omega0inv %*% G %*% H0
  
  # compute optimal bandwidth
  nu2 <- (2 * kernConst$mu1 + kernConst$mu2) * (k - p) * sum(diag(Sigma0 %*% W))
  nu3 <- kernConst$g_q^2 * sum(diag(t(Omega_q) %*% t(H0) %*% W %*% H0 %*% Omega_q %*% P0))
  c0 <- if(nu2 * nu3 > 0) 2 * kernConst$q else -1
  rval <- (c0 * nu3/nu2 * n)^(1/(1 + 2 * kernConst$q))

  return(rval)
}



