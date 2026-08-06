#' ARMA copula processes
#'
#' Class of objects for ARMA copula processes.
#'
#' @slot name name of ARMA copula process.
#' @slot modelspec vector containing number of AR and MA parameters.
#' @slot pars list consisting of vector of AR parameters named `ar`
#' and vector of MA parameters named `ma`.
#'
#' @export
#'
setClass("armacopula", contains = "tscopula", slots = list(
  name = "character",
  modelspec = "numeric",
  pars = "list"
))

#' Constructor function for ARMA copula process
#' @param pars list consisting of vector of AR parameters named `ar`
#' and vector of MA parameters named `ma`.
#'
#' @return An object of class \linkS4class{armacopula}.
#' @export
#'
#' @examples
#' armacopula(list(ar = 0.5, ma = 0.4))
armacopula <- function(pars = list(ar = 0, ma = 0)) {
  if ("ar" %in% names(pars)) {
    arpars <- pars$ar
    if (is.null(arpars)) {
      stop("No NULL values for parameters; omit from list instead")
    }
    if (non_stat(arpars)) {
      stop("Non-stationary AR model")
    }
    p <- length(arpars)
    names(pars$ar) <- paste("ar", 1:p, sep = "")
  }
  else {
    p <- 0
  }
  if ("ma" %in% names(pars)) {
    mapars <- pars$ma
    if (is.null(mapars)) {
      stop("No NULL values for parameters; omit from list instead")
    }
    if (non_invert(mapars)) {
      stop("Non-invertible MA model")
    }
    q <- length(mapars)
    names(pars$ma) <- paste("ma", 1:q, sep = "")
  }
  else {
    q <- 0
  }
  if ((p == 0) & (q == 0)) {
    stop("Specify named ar and/or ma parameters")
  }
  new("armacopula",
      name = paste("ARMA(", p, ",", q, ")", sep = ""),
      modelspec = c(p = p, q = q),
      pars = pars
  )
}

#' @describeIn armacopula Coef method for ARMA copula class
#'
#' @param object an object of the class.
#'
#' @export
setMethod("coef", "armacopula", function(object) {
  p <- object@modelspec[1]
  q <- object@modelspec[2]
  if (p > 0) {
    arpars <- object@pars$ar
  } else {
    arpars <- NULL
  }
  if (q > 0) {
    mapars <- object@pars$ma
  } else {
    mapars <- NULL
  }
  c(arpars, mapars)
})

#' @describeIn armacopula Show method for ARMA copula process
#'
#' @param object an object of the class.
#'
#' @export
#'
setMethod("show", c(object = "armacopula"), function(object) {
  cat("object class: ", is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: \n")
  print(coef(object))
})

#' Check for causality of ARMA process
#'
#' @param ar vector of autoregressive parameters
#'
#' @return A logical variable stating whether ARMA process is causal.
#' @export
#'
non_stat <- function(ar) {
  status <- FALSE
  if (sum(ar^2) > 0) {
    roots <- polyroot(c(1, -ar))
    if (min(abs(roots)) <= 1) {
      status <- TRUE
    }
  }
  status
}

#' Check for invertibility of ARMA process
#'
#' @param ma vector of moving average parameters.
#'
#' @return A logical variable stating whether ARMA process is invertible.
#' @export
#'
non_invert <- function(ma) {
  status <- FALSE
  if (sum(ma^2) > 0) {
    roots <- polyroot(c(1, ma))
    if (min(abs(roots)) <= 1) {
      status <- TRUE
    }
  }
  status
}

#' @describeIn armacopula Simulation method for armacopula class
#'
#' @param object an object of the class.
#' @param n length of realization.
#'
#' @export
#'
#' @examples
#' sim(armacopula(list(ar = c(0.5, 0.4), ma = -0.8)), n = 1000)
setMethod("sim", c(object = "armacopula"), function(object, n = 1000) {
  discard <- sum(object@modelspec)
  pnorm(arima.sim(
    model = object@pars,
    n = n,
    n.start = discard + 1,
    sd = sigmastarma(object)
  ))
})

#' Standard deviation of innovations for armacopula
#'
#' Uses the function \code{\link[ltsa]{tacvfARMA}} in the ltsa library.
#'
#' @param x an object of class \linkS4class{armacopula}.
#'
#' @return The standard deviation of the standardized ARMA innovation distribution.
#' @export
#'
#' @examples
#' sigmastarma(armacopula(list(ar = c(0.5, 0.4), ma = -0.8)))
sigmastarma <- function(x) {
  ar <- x@pars$ar
  ma <- x@pars$ma
  if (length(ar) == 0) {
    ar <- 0
  }
  if (length(ma) == 0) {
    ma <- 0
  }
  1 / sqrt(ltsa::tacvfARMA(phi = ar, theta = -ma, maxLag = 0, sigma2 = 1))
}

#' Objective function for ARMA copula process
#'
#' @param theta vector of parameters of ARMA process
#' @param modelspec vector containing model order (p,q)
#' @param u vector of data
#'
#' @return Value of objective function at parameters.
#' @keywords internal
#'
armacopula_objective <- function(theta, modelspec, u) {
  xdata <- qnorm(u)
  p <- modelspec[1]
  ar <- 0
  q <- modelspec[2]
  ma <- 0
  if (p > 0) {
    ar <- theta[1:p]
  }
  if (q > 0) {
    ma <- theta[(p + 1):(p + q)]
  }
  if (non_stat(ar) | non_invert(ma)) {
    output <- NA
  } else {
    sp <- starmaStateSpace(ar, ma, c(p, q))
    ans <- FKF::fkf(
      a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt, Zt = sp$Zt,
      HHt = sp$HHt, GGt = sp$GGt, yt = rbind(xdata)
    )
    output <- -ans$logLik + sum(log(dnorm(xdata)))
  }
  return(output)
}

#' State space representation for standardized ARMA model
#'
#' @param ar vector of ar parameters
#' @param ma vector of ma parameters
#' @param order vector giving order (p,q)
#'
#' @return State space representation of ARMA process.
#' @keywords internal
#'
starmaStateSpace <- function(ar, ma, order) {
  p <- order[1]
  q <- order[2]
  m <- max(p, q + 1)
  allar <- rep(0, m)
  allma <- rep(0, m)
  if (p > 0) {
    allar[1:p] <- ar
  }
  if (q > 0) {
    allma[1:q] <- ma
  }
  Tt <- matrix(allar)
  Zt <- matrix(1)
  if (m > 1) {
    block1 <- diag(m - 1)
    block2 <- matrix(0, ncol = m - 1)
    rmat <- rbind(block1, block2)
    Tt <- cbind(Tt, rmat)
    Zt <- cbind(Zt, matrix(0, ncol = m - 1))
  }
  ct <- matrix(0)
  dt <- matrix(0, nrow = m)
  GGt <- matrix(0)
  sigma2 <- ltsa::tacvfARMA(phi = ar, theta = -ma, maxLag = 0, sigma2 = 1)
  Hcontent <- 1
  if (m > 1) {
    Hcontent <- c(1, allma[1:(m - 1)])
  }
  H <- matrix(Hcontent, nrow = m) / sqrt(sigma2)
  HHt <- H %*% t(H)
  ## initialization
  a0 <- rep(0, m)
  P0 <- matrix(solve(diag(1, m^2) - kronecker(Tt, Tt)) %*% as.vector(HHt),
    nrow = m,
    ncol = m
  )
  return(list(a0 = a0, P0 = P0, ct = ct, dt = dt, Zt = Zt, Tt = Tt, GGt = GGt, HHt = HHt))
}

#' Kalman filter for ARMA copula model
#'
#' @param x an object of class \linkS4class{armacopula}.
#' @param y a vector of data.
#'
#' @return A matrix or multivariate time series with columns consisting of
#' conditional mean, standard deviation and residuals.
#' @export
#'
#' @examples
#' data <- sim(armacopula(list(ar = c(0.5, 0.4), ma = -0.8)), n = 1000)
#' kfilter(armacopula(list(ar = c(0.5, 0.4), ma = -0.8)), data)
kfilter <- function(x, y) {
  n <- length(y)
  ar <- x@pars$ar
  p <- x@modelspec[1]
  if (p == 0) {
    ar <- 0
  }
  ma <- x@pars$ma
  q <- x@modelspec[2]
  if (q == 0) {
    ma <- 0
  }
  sp <- starmaStateSpace(ar, ma, c(p, q))
  ans <- FKF::fkf(
    a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt, Zt = sp$Zt,
    HHt = sp$HHt, GGt = sp$GGt, yt = rbind(as.numeric(qnorm(y)))
  )
  mu_t <- ans$at[1, 1:n]
  sigma_t <- sqrt(ans$Ft[1, 1, 1:n])
  resid <- ans$vt[1, 1:n]
  if (inherits(x, "zoo")) {
    attributes(resid) <- attributes(x)
    attributes(mu_t) <- attributes(x)
    attributes(sigma_t) <- attributes(x)
  }
  fseries <- cbind(mu_t, sigma_t, resid)
  dimnames(fseries) <- list(NULL, c("mu_t", "sigma_t", "resid"))
  fseries
}

#' Residual function for armacopula object
#'
#' @param object a fitted armacopula object.
#' @param data the data to which copula is fitted.
#' @param trace extract trace instead of residuals.
#'
#' @return vector of model residuals
#' @keywords internal
#'
resid_armacopula <- function(object, data = NA, trace = FALSE){
  series <- kfilter(object, data)
  if (trace)
    output <- series[, "mu_t"]
  else
    output <- series[, "resid"]/sigmastarma(object)
  output
}

#' @describeIn armacopula Calculate Kendall's tau values for armacopula model
#'
#' @param object an object of the class.
#' @param lagmax maximum value of lag.
#'
#' @export
#'
#' @examples
#' mod <- armacopula(list(ar = 0.95, ma = -0.85))
#' kendall(mod)
setMethod("kendall", c(object = "armacopula"), function(object, lagmax = 20){
  ar <- 0
  ma <- 0
  if (object@modelspec[1] > 0)
    ar <- object@pars$ar
  if (object@modelspec[2] > 0)
    ma <- object@pars$ma
  pacf <- ARMAacf(ar = ar, ma = ma, lag.max = lagmax, pacf = TRUE)
  tau <- (2/pi)*asin(pacf)
  tau
}
)

#' Generalized lagging for fitted armacopula objects
#'
#' @param copula an armacopula object.
#' @param data the data to which copula is fitted.
#' @param lagmax the maximum lag value.
#' @param glagplot logical value indicating generalized lag plot.
#'
#' @return If \code{glagplot} is \code{TRUE} a list of generalized lagged datasets
#' of maximum length 9 is returned to facilitate a generalized lagplot.
#' If \code{glagplot} is \code{FALSE} a vector of length \code{lagmax} containing
#' the Kendall rank correlations for the generalized lagged datasets is returned.
#'
#' @keywords internal
glag_for_armacopula <- function(copula, data, lagmax, glagplot = FALSE) {
  n <- length(data)
  k <- lagmax
  data <- cbind(as.numeric(data[1:(n - 1)]), as.numeric(data[2:n]))
  if (glagplot){
    k <- min(k, 9)
    output <- vector(mode = "list", length = k)
    output[[1]] <- data
  }
  else{
  output <- rep(NA, k)
  output[1] <- cor(data, method = "kendall")[1, 2]
  }
  ar <- 0
  ma <- 0
  if (copula@modelspec[1] > 0)
    ar <- copula@pars$ar
  if (copula@modelspec[2] > 0)
    ma <- copula@pars$ma
  pacf <- ARMAacf(ar = ar, ma = ma, lag.max = lagmax, pacf = TRUE)
  if (k >1){
    for (i in 1:(k - 1)) {
    n <- dim(data)[1]
    model <- rvinecopulib::bicop_dist(family = "gauss", parameters = pacf[i])
    data <-
      cbind(rvinecopulib::hbicop(data[(1:(n - 1)), ], model, cond_var = 2),
            rvinecopulib::hbicop(data[(2:n), ], model, cond_var = 1))
    if (glagplot)
      output[[i+1]] <- data
    else
      output[i+1] <- cor(data, method = "kendall")[1, 2]
    }
  }
output
}

#' @describeIn armacopula Prediction method for armacopula class
#'
#' @param object an object of the class.
#' @param data vector of past data values.
#' @param x vector of arguments of prediction function.
#' @param type type of prediction function ("df" for density, "qf" for quantile function
#' or "dens" for density).
#'
#' @export
#'
setMethod("predict", c(object = "armacopula"), function(object, data, x, type = "df") {
  ar <- object@pars$ar
  p <- object@modelspec[1]
  if (p == 0) {
    ar <- 0
  }
  ma <- object@pars$ma
  q <- object@modelspec[2]
  if (q == 0) {
    ma <- 0
  }
  sp <- starmaStateSpace(ar, ma, c(p, q))
  ans <- FKF::fkf(
    a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt, Zt = sp$Zt,
    HHt = sp$HHt, GGt = sp$GGt, yt = rbind(as.numeric(qnorm(data)))
  )
  n <- length(data)
  mu_t <- ans$at[1, (n+1)]
  sigma_t <- sqrt(ans$Pt[1, 1, (n+1)])
  switch(type,
         "df" = pnorm(qnorm(x), mean = mu_t, sd = sigma_t),
         "qf" = pnorm(qnorm(x, mean = mu_t, sd = sigma_t)),
         "dens" = dnorm(qnorm(x), mean = mu_t, sd = sigma_t)/dnorm(qnorm(x)))
})

#' Transform an armacopula into a dvinecopula or dvinecopula2 object
#'
#' @param object an object of class \linkS4class{armacopula}.
#'
#' @return An object of class \linkS4class{dvinecopula} (for AR copulas)
#' or class \linkS4class{dvinecopula2} (for MA or ARMA copulas).
#' @export
#'
#' @examples
#' arma2dvine(armacopula(list(ar = 0.5, ma = 0.4)))
arma2dvine <- function(object){
  if (!(is(object, "armacopula")))
    stop("Not armacopula object")
  ord <- object@modelspec
  p <- ord[1]
  q <- ord[2]
  if (q == 0)
    tscop <- dvinecopula(family = "gauss", pars = coef(object))
  else
    tscop <- dvinecopula2(family = "gauss", pars = object@pars)
  tscop
}
