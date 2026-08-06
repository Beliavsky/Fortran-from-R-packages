#' D-vine copula processes of type 2
#'
#' Class of objects for d-vine copula processes. See \link{dvinecopula2} for more details.
#'
#' @slot name name of the d-vine copula process.
#' @slot modelspec list containing the family, rotation, and name of KPACF
#' @slot pars list comprising of the parameters.
#'
#' @export
#'
setClass("dvinecopula2", contains = "tscopula", slots = list(
  name = "character",
  modelspec = "list",
  pars = "list"
))

#' Constructor function for dvinecopula2 process
#'
#' This function sets up a stationary d-vine process of finite or infinite order based on a single
#' copula family from a subset of those that can be implemented using
#' \code{\link[rvinecopulib]{bicop_dist}} in the \code{rvinecopulib} package.
#'
#' The copula family may be any one-parameter family or the t copula family. The basic copula from
#' which the sequence is built may be rotated through 180 degrees using the \code{rotation} argument; the default
#' is no rotation (0 degrees).
#'
#' The copulas are parameterized using the Kendall partial autocorrelation function (kpacf) specified
#' by the \code{kpacf} argument. The default choice is the kpacf of a standard ARMA process which is
#' implemented in the function \code{\link{kpacf_arma}}. The parameters
#' of the kpacf should be set as a list using the \code{pars} argument; the required parameters should usually
#' be clear from the documentation of the chosen kpacf function and must be correctly named.
#'
#' If the kpacf takes a negative value at any lag and the standard copula is unable to model a
#' negative dependency (e.g. Clayton, Gumbel, Joe and their 180 degree rotations) then one of four
#' different treatments may be specified using the \code{negtau} parameter: "gauss" substitutes a
#' Gaussian copula at that lag; "frank" substitutes a Frank copula; "right" and "left" rotate the copula
#' through 90 degrees in a clockwise or anto-clockwise direction respectively.
#'
#' The \code{maxlag} parameter specifies the maximum lag of the process; a finite number gives a finite-order
#' stationary d-vine process, but the parameter may also be set to \code{Inf} for an infinite-order process.
#'
#' If the t copula is chosen by setting \code{family} equal to "t", the list of
#' parameters needs to be augmented with a component named "df" which is
#' the degrees of freedom. In this case it makes sense to set \code{maxlag} to be a finite number to avoid models
#' with tail dependencies at arbitrary lags which are not ergodic. The class \linkS4class{dvinecopula3}
#' is more suitable for working with t copulas with different degrees of freedom at different lags.
#'
#' @param family family name
#' @param rotation a scalar specifying the rotation (default is 0)
#' @param kpacf a character string giving the name of the Kendall pacf
#' @param pars a list containing the parameters of the model
#' @param maxlag a scalar specifying the maximum lag
#' @param negtau a character string specifying the treatment of negative Kendall's tau values
#'
#' @return An object of class \linkS4class{dvinecopula2}.
#' @export
#'
#' @examples
#' dvinecopula2(family = "joe", kpacf = "kpacf_arma",
#' pars = list(ar = 0.95, ma = -0.85), maxlag = 30)
dvinecopula2 <- function(family = "gauss",
                         rotation = 0,
                         kpacf = "kpacf_arma",
                         pars = list(ar = 0.1, ma = 0.1),
                         maxlag = Inf,
                         negtau = "none") {
  if (!(is(family, "character")))
    stop("copula family must be specified by name")
  if (is.null(names(pars)))
    stop("parameters should be named (p1 and p2 for exp/power)")
  fam <- tolower(family)
  if (fam %in% c("gauss", "frank", "t"))
    negtau <- "none"
 modelspec <- list(family = fam,
                   rotation = rotation,
                   kpacf = kpacf,
                   maxlag = maxlag,
                   npar = length(unlist(pars)),
                   negtau = negtau)
  new("dvinecopula2",
      name = paste("type2-d-vine"),
      modelspec = modelspec,
      pars = pars
  )
}

#' KPACF of ARMA process
#'
#' @param k number of lags.
#' @param theta list with components ar and ma specifying the ARMA parameters.
#'
#' @return A vector of Kendall partial autocorrelations of length \code{k}.
#' @export
#'
kpacf_arma <- function(k, theta){
  if (is.list(theta))
    theta <- unlist(theta)
  ar <- numeric()
  ma <- numeric()
  nm <- substring(names(theta), 1, 2)
  if ("ar" %in% nm)
    ar <- theta[nm == "ar"]
  if ("ma" %in% nm)
    ma <- theta[nm == "ma"]
  if ((non_stat(ar)) | (non_invert(ma)))
    return(rep(NA, k))
  pacf <- ARMAacf(ar = ar, ma = ma, lag.max = k, pacf = TRUE)
  (2/pi)*asin(pacf)
}

#' KPACF of quarterly seasonal ARMA process
#'
#' @param k number of lags.
#' @param theta list with components ar, ma, sar and sma specifying the ARMA and seasonal ARMA parameters.
#'
#' @return A vector of Kendall partial autocorrelations of length \code{k}.
#' @export
kpacf_sarma4 <- function (k, theta)
{
  if (is.list(theta))
    theta <- unlist(theta)
  nm <- substring(names(theta), 1, 2)
  D <- 4
  ar <- numeric()
  ma <- numeric()
  sar <- numeric()
  sma <- numeric()
  if ("ar" %in% nm)
    ar <- theta[nm == "ar"]
  if ("ma" %in% nm)
    ma <- theta[nm == "ma"]
  nm <- substring(names(theta), 1, 3)
  if ("sar" %in% nm){
    sar <- theta[nm == "sar"]
    sar <- as.vector(sapply(sar, function(x, d){c(rep(0, d-1), x)}, d=D))
    ar <- -coefficients(polynom::polynomial(c(1, -sar)) * polynom::polynomial(c(1, -ar)))[-1]
    if (length(ar) == 0) ar <-0
  }
  if ("sma" %in% nm){
    sma <- theta[nm == "sma"]
    sma <- as.vector(sapply(sma, function(x, d){c(rep(0, d-1), x)}, d=D))
    ma <- coefficients(polynom::polynomial(c(1, sma)) * polynom::polynomial(c(1, ma)))[-1]
    if (length(ma) == 0) ma <- 0
  }
  if ((non_stat(ar)) | (non_invert(ma)))
    return(rep(NA, k))
  pacf <- ARMAacf(ar = ar, ma = ma, lag.max = k, pacf = TRUE)
  (2/pi) * asin(pacf)
}

#' KPACF of monthly seasonal ARMA process
#'
#' @param k number of lags.
#' @param theta list with components ar, ma, sar and sma specifying the ARMA and seasonal ARMA parameters.
#'
#' @return A vector of Kendall partial autocorrelations of length \code{k}.
#' @export
kpacf_sarma12 <- function (k, theta)
{
  if (is.list(theta))
    theta <- unlist(theta)
  nm <- substring(names(theta), 1, 2)
  D <- 12
  ar <- numeric()
  ma <- numeric()
  sar <- numeric()
  sma <- numeric()
  if ("ar" %in% nm)
    ar <- theta[nm == "ar"]
  if ("ma" %in% nm)
    ma <- theta[nm == "ma"]
  nm <- substring(names(theta), 1, 3)
  if ("sar" %in% nm){
    sar <- theta[nm == "sar"]
    sar <- as.vector(sapply(sar, function(x, d){c(rep(0, d-1), x)}, d=D))
    ar <- -coefficients(polynom::polynomial(c(1, -sar)) * polynom::polynomial(c(1, -ar)))[-1]
    if (length(ar) == 0) ar <-0
  }
  if ("sma" %in% nm){
    sma <- theta[nm == "sma"]
    sma <- as.vector(sapply(sma, function(x, d){c(rep(0, d-1), x)}, d=D))
    ma <- coefficients(polynom::polynomial(c(1, sma)) * polynom::polynomial(c(1, ma)))[-1]
    if (length(ma) == 0) ma <- 0
  }
  if ((non_stat(ar)) | (non_invert(ma)))
    return(rep(NA, k))
  pacf <- ARMAacf(ar = ar, ma = ma, lag.max = k, pacf = TRUE)
  (2/pi) * asin(pacf)
}

#' Compute partial autocorrelations from autocorrelations
#'
#' @param rho vector of autocorrelation values (excluding 1).
#'
#' @return A vector of partial autocorrelation values with same length as \code{rho}.
#' @export
#'
#' @examples
#' rho <- ARMAacf(ar = -0.9, ma = 0.8, lag.max = 50)[-1]
#' alpha <- acf2pacf(rho)
acf2pacf <- function(rho) {
  L <- length(rho)
  pi <- numeric(L)
  pi[1] <- rho[1]
  phik <- pi
  vk <- (1 - pi[1] ^ 2)
  if (L > 1) {
    for (k in 2:L) {
      a <- sum(c(1, -phik[1:(k - 1)]) * rev(rho[1:k])) / vk
      phik <- c(phik[1:(k - 1)] - a * rev(phik[1:(k - 1)]), a)
      vk <- vk * (1 - a ^ 2)
      pi[k] <- a
    }
  }
  pi
}

#' Compute autocorrelations from partial autocorrelations
#'
#' @param alpha vector of partial autocorrelation values.
#'
#' @return A vector of autocorrelation values with same length as \code{alpha}.
#' @export
#'
#' @examples
#' alpha <- ARMAacf(ar = -0.9, ma = 0.8, lag.max = 50, pacf = TRUE)
#' rho <- pacf2acf(alpha)
pacf2acf <- function(alpha) {
  L <- length(alpha)
  phik <- numeric(L)
  phik[1] <- alpha[1]
  if (L > 1)
    for (k in 2:L)
      phik[1:k] = c(phik[1:(k - 1)] - alpha[k] * rev(phik[1:(k - 1)]), alpha[k])
  rho <- stats::ARMAacf(ar = phik, lag.max = length(alpha))
  as.numeric(rho[-1])
}

#' Compute autoregressive coefficients from partial autocorrelations
#'
#' @param alpha vector of partial autocorrelation values.
#'
#' @return A vector of autoregressive coefficients with same length as \code{alpha}.
#' @export
#'
#' @examples
#' alpha <- ARMAacf(ar = -0.9, ma = 0.8, lag.max = 50, pacf = TRUE)
#' phi <- pacf2ar(alpha)
pacf2ar <- function(alpha) {
  L <- length(alpha)
  phik <- numeric(L)
  phik[1] <- alpha[1]
  if (L > 1)
    for (k in 2:L)
      phik[1:k] = c(phik[1:(k - 1)] - alpha[k] * rev(phik[1:(k - 1)]), alpha[k])
  phik
}

#' KPACF of ARFIMA process
#'
#' @param k number of lags.
#' @param theta list with components ar, ma and d specifying the ARFIMA parameters
#'
#' @return A vector of Kendall partial autocorrelations of length \code{k}.
#' @export
#'
kpacf_arfima <- function(k, theta){
  if (is.list(theta))
    theta <- unlist(theta)
  ar <- numeric()
  ma <- numeric()
  nm <- substring(names(theta), 1, 2)
  if ("ar" %in% nm)
    ar <- theta[nm == "ar"]
  if ("ma" %in% nm)
    ma <- theta[nm == "ma"]
  if ("d" %in% nm)
    d <- theta[nm == "d"]
  if ((non_stat(ar)) | (non_invert(ma)) | (d <= -0.5) | (d >= 0.5))
    return(rep(NA, k))
  acvf <- arfima::tacvfARFIMA(phi = ar, theta = -ma, dfrac = d, maxlag = k)
  acf <- acvf[-1]/acvf[1]
  pacf <- acf2pacf(acf)
  (2/pi)*asin(pacf)
}

#' KPACF of fractional Brownian noise
#'
#' @param k number of lags
#' @param theta parameter of process
#'
#' @return A vector of Kendall partial autocorrelations of length \code{k}.
#' @export
#'
kpacf_fbn <- function(k, theta){
  if (is.list(theta))
    theta <- unlist(theta)
  if ((theta <= 0) | (theta >= 1))
    return(rep(NA, k))
  acf <- (((1:k) + 1)^{2 * theta[1]} + abs((1:k) - 1)^{2 * theta[1]} - 2 * (1:k)^{2 * theta[1]})/2
  pacf <- acf2pacf(acf)
  return((2/pi)*asin(pacf))
}

#' @describeIn dvinecopula2 Coef Method for dvinecopula2 class
#'
#' @param object an object of the class.
#'
#' @export
#'
setMethod("coef", c(object = "dvinecopula2"), function(object) {
  if (length(object@pars) == 1) {
    return(object@pars[[1]])
  } else {
    nms <- unlist(lapply(object@pars, names), use.names = FALSE)
    vals <- unlist(object@pars, use.names = FALSE)
    names(vals) <- nms
    return(vals)
  }
})

#' @describeIn dvinecopula2 Show method for dvinecopula2 class
#'
#' @param object an object of the class.
#'
#' @export
#'
setMethod("show", c(object = "dvinecopula2"), function(object) {
  cat("object class: ", is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  famname <- object@modelspec$family
  if (object@modelspec$rotation !=0)
    famname <- paste(famname, "with rotation", object@modelspec$rotation)
  cat("copula family: ", famname, "\n", sep = "")
  if (object@modelspec$negtau != "none")
    cat("negative tau treatment: ", object@modelspec$negtau, "\n", sep = "")
  kpacf  <- object@modelspec$kpacf
  if (object@modelspec$maxlag != Inf)
    kpacf <- paste(kpacf, "with max lag", object@modelspec$maxlag)
  cat("KPACF: ", kpacf,"\n", sep = "")
  cat("parameters: \n")
  print(coef(object))
})

#' Objective function for dvinecopula2 process
#'
#' @param theta parameters of kpacf
#' @param modelspec list specifying model
#' @param u data
#'
#' @return Value of objective function at parameters.
#'
#' @keywords internal
#'
dvinecopula2_objective <- function(theta, modelspec, u) {
  n <- length(u)
  kpacf <- eval(parse(text = modelspec$kpacf))
  tauvals <- kpacf((n-1), theta)
  if (is.na(sum(tauvals)))
    return(NA)
  k <- 1
  largetau <- (abs(tauvals) > .Machine$double.eps)
  if (sum(largetau) > 0)
    k <- max((1:(n-1))[largetau])
  k <- min(k, modelspec$maxlag)
  pc_list <- vector("list", k)
  for (i in 1:k) {
    fam <- tolower(modelspec$family)
    rot <- modelspec$rotation
    if (tauvals[i] < 0){
      if (modelspec$negtau == "left")
        rot <- rot + 90
      if (modelspec$negtau == "right")
        rot <- (rot + 270) %% 360
      if (modelspec$negtau %in% c("gauss","frank")){
        fam <- modelspec$negtau
        rot <- 0
      }
      if ((modelspec$negtau == "none") & (fam %in% c("gumbel", "joe", "clayton")))
        return(NA)
    }
    coppars <- ktau_to_par(
      family = fam,
      tau = tauvals[i]
    )
    if (is.na(coppars))
      return(NA)
    if (fam == "t")
      coppars <- c(coppars, theta["df"])
    pc_list[[i]] <- tryCatch(rvinecopulib::bicop_dist(
      family = fam,
      rotation = rot,
      parameters = coppars
    ),
    error = function(e) {
      return(NA)
    }
    )
    if (is.na(pc_list[[i]][[1]])) {
      return(NA)
    }
  }
  v <- cbind(u[1:(n - 1)], u[2:n])
  LL <- 0
  for (j in 1:k) {
    LL <- LL + sum(log(rvinecopulib::dbicop(u = v, family = pc_list[[j]])))
    if (j == k) {
      return(-LL)
    }
    n <- dim(v)[1]
    v <- cbind(
      rvinecopulib::hbicop(v[(1:(n - 1)), ], cond_var = 2, family = pc_list[[j]]),
      rvinecopulib::hbicop(v[(2:n), ], cond_var = 1, family = pc_list[[j]])
    )
  }
}

#' Transform Kendall's tau values to copula parameters
#'
#' @param family name of copula family
#' @param tau value of Kendall's tau
#'
#' @return Scalar value of parameter giving tau for a copula family.
#' @keywords internal
ktau_to_par <- function(family, tau){
  if (family == "t")
    family <- "gauss"
  if (family == "bb1")
    family <- "clayton"
  rvinecopulib::ktau_to_par(family, tau)
}

#' @describeIn dvinecopula2 Simulation method for dvinecopula2 class
#'
#' @param object an object of the class.
#' @param n length of realization.
#'
#' @export
#'
setMethod("sim", c(object = "dvinecopula2"), function(object, n = 1000) {
  pc_list <- mklist_dvine2(object, n-1, truncate = TRUE, tol = 1/3)
  simdvine(pc_list, n, innov = NA, start = NA)
})

#' Make list of pair copulas for dvinecopula2 object
#'
#' @param x an object of class dvinecopula2
#' @param maxlag maximum possible lag to consider
#' @param truncate logical variable stating whether to truncate copulas with negligible dependence
#' @param tol tolerance for truncation
#'
#' @return a list of pair copulas
#' @keywords internal
#'
mklist_dvine2 <- function(x, maxlag, truncate, tol = 1){
  k <- min(maxlag, x@modelspec$maxlag)
  kpacf <- eval(parse(text = x@modelspec$kpacf))
  tauvals <- kpacf(k, x@pars)
  if (truncate)
    k <- max(c(1, which(abs(tauvals) > .Machine$double.eps^tol)))
  pc_list <- vector("list", k)
  for (i in 1:k) {
    fam <- tolower(x@modelspec$family)
    rot <- x@modelspec$rotation
    if (tauvals[i] < 0){
      if (x@modelspec$negtau == "left")
        rot <- rot + 90
      if (x@modelspec$negtau == "right")
        rot <- (rot + 270) %% 360
      if (x@modelspec$negtau %in% c("gauss","frank")){
        fam <- x@modelspec$negtau
        rot <- 0
      }
    }
    coppars <- ktau_to_par(
      family = fam,
      tau = tauvals[i]
    )
    if (fam == "t")
      coppars <- c(coppars, x@pars$df)
    pc_list[[i]] <- rvinecopulib::bicop_dist(
      family = fam,
      rotation = rot,
      parameters = coppars)
  }
  pc_list
}

#' @describeIn dvinecopula2 Prediction method for dvinecopula2 class
#'
#' @param object an object of the class.
#' @param data vector of past data values.
#' @param x vector of arguments of prediction function.
#' @param type type of prediction function ("df" for density, "qf" for quantile function
#' or "dens" for density).
#'
#' @export
#'
setMethod("predict", c(object = "dvinecopula2"), function(object, data, x, type = "df") {
  pc_list <- mklist_dvine2(object, length(data)-1, truncate = TRUE, tol = 1/3)
  switch(type,
         "df" = Rblatt(pc_list, data, x),
         "qf" = IRblatt(pc_list, data, x),
         "dens" = Rblattdens(pc_list, data, x))

})

#' Residual function for dvinecopula2 object
#'
#' @param object a fitted dvinecopula2 object.
#' @param data the data to which copula is fitted.
#' @param trace extract trace instead of residuals.
#'
#' @return vector of model residuals
#' @keywords internal
#'
resid_dvinecopula2 <- function(object, data = NA, trace = FALSE){
  n <- length(data)
  pc_list <- mklist_dvine2(object, n-1, truncate = TRUE, tol = 1/3)
  k <- length(pc_list)
  for (i in 1:k)
    if (pc_list[[i]]$rotation %in% c(90,270))
      pc_list[[i]]$rotation <- 360 - pc_list[[i]]$rotation
  if (trace)
    target <- rep(0.5, n)
  else
    target <- data
  res <- rep(NA, n)
  res[1] <- target[1]
  if (k >1){
    for (i in 2:k){
      pcs <- lapply(1:(i-1), function(j) {
        replicate(i - j, pc_list[[j]], simplify = FALSE)
      })
      vc_short <- rvinecopulib::vinecop_dist(pcs, rvinecopulib::dvine_structure(i:1))
      vals <- c(data[1:(i-1)], target[i])
      res[i] <- rvinecopulib::rosenblatt(t(vals), vc_short)[i]
    }
  }
  pcs <- lapply(1:k, function(j) {
    replicate(k + 1 - j, pc_list[[j]], simplify = FALSE)
  })
  vc_short <- rvinecopulib::vinecop_dist(pcs, rvinecopulib::dvine_structure((k + 1):1))
  for (i in ((k+1):n)){
    vals <- c(data[(i-k):(i-1)], target[i])
    res[i] <- rvinecopulib::rosenblatt(t(vals), vc_short)[k+1]
  }
  qnorm(res)
}

#' @describeIn dvinecopula2 Calculate Kendall's tau values for pair copulas in type 2 d-vine copula
#'
#' @param object an object of the class.
#' @param lagmax maximum value of lag.
#'
#' @export
#'
#' @examples
#' copmod <- dvinecopula2(family = "joe", kpacf = "kpacf_arma",
#' pars = list(ar = 0.95, ma = -0.85), maxlag = 30)
#' kendall(copmod)
setMethod("kendall", c(object = "dvinecopula2"), function(object, lagmax = 20) {
  kpacf <- eval(parse(text = object@modelspec$kpacf))
  kpacf(lagmax, object@pars)
}
)

#' Generalized lagging for fitted dvinecopula2 objects
#'
#' @param copula a dvinecopula2 object
#' @param data the data to which copula is fitted
#' @param lagmax the maximum lag value.
#' @param glagplot logical value indicating generalized lag plot.
#'
#' @return If \code{glagplot} is \code{TRUE} a list of generalized lagged datasets
#' of maximum length 9 is returned to facilitate a generalized lagplot.
#' If \code{glagplot} is \code{FALSE} a vector of length \code{lagmax} containing
#' the Kendall rank correlations for the generalized lagged datasets is returned.
#' @keywords internal
glag_for_dvinecopula2 <- function(copula, data, lagmax, glagplot = FALSE) {
  if (glagplot)
    lagmax <- min(lagmax, 9)
  pc_list <- mklist_dvine2(copula, lagmax, truncate = FALSE)
  k <- length(pc_list)
  n <- length(data)
  data <- cbind(as.numeric(data[1:(n - 1)]), as.numeric(data[2:n]))
  if (glagplot){
    output <- vector(mode = "list", length = k)
    output[[1]] <- data
  }
  else{
    output <- rep(NA, k)
    output[1] <- cor(data, method = "kendall")[1, 2]
  }
  if (k >1){
    for (i in 1:(k - 1)) {
      n <- dim(data)[1]
      data <-
        cbind(rvinecopulib::hbicop(data[(1:(n - 1)), ], pc_list[[i]], cond_var = 2),
              rvinecopulib::hbicop(data[(2:n), ], pc_list[[i]], cond_var = 1))
      if (glagplot)
        output[[i+1]] <- data
      else
        output[i+1] <- cor(data, method = "kendall")[1, 2]
    }
  }
  output
}

