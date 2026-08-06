#' Time series copula processes with v-transforms
#'
#' Class of objects for v-transformed time series copula processes.
#'
#' @slot Vcopula object of class \linkS4class{tscopulaU}.
#' @slot Vtransform object of class \linkS4class{Vtransform}.
#' @slot Wcopula object of class \linkS4class{tscopula}.
#'
#' @export
#'
setClass("vtscopula",
  contains = "tscopula",
  slots = list(
    Vcopula = "tscopulaU",
    Vtransform = "Vtransform",
    Wcopula = "tscopula"
  )
)

#' Constructor function for vtscopula object
#'
#' @param tscopulaU an object of class
#' \linkS4class{armacopula}, \linkS4class{dvinecopula} or \linkS4class{dvinecopula2}.
#' @param Vtransform an object of class \linkS4class{Vtransform}.
#' @param Wcopula an object of class \linkS4class{tscopula}.
#'
#' @return An object of class \linkS4class{vtscopula}.
#' @export
#'
#' @examples
#' copobject <- armacopula(pars = list(ar = 0.6, ma = 0.2))
#' vtscopula(copobject, Vtransform = V2p())
vtscopula <- function(tscopulaU,
                      Vtransform = Vlinear(),
                      Wcopula = swncopula()) {
  new("vtscopula",
    Vcopula = tscopulaU,
    Vtransform = Vtransform,
    Wcopula = Wcopula
  )
}

#' @describeIn vtscopula Show method for vtscopula objects
#'
#' @param object an object of the class.
#'
#' @export
#'
setMethod("show", "vtscopula", function(object) {
  cat("object class: ", is(object)[[1]], "\n", sep = "")
  cat("____________ \n")
  cat("Base copula: \n")
  show(object@Vcopula)
  cat("____________ \n")
  cat("V-transform: \n")
  show(object@Vtransform)
  if (!(is(object@Wcopula, "swncopula"))){
    cat("________ \n")
    cat("Wcopula: \n")
    show(object@Wcopula)
  }
})

#' @describeIn vtscopula Coef method for vtscopula class
#'
#' @param object an object of the class.
#'
#' @export
#'
setMethod("coef", "vtscopula", function(object) {
  c(coef(object@Vcopula), coef(object@Vtransform), coef(object@Wcopula))
})

#' @describeIn vtscopula Prediction method for vtscopula class
#'
#' @param object an object of the class.
#' @param data vector of past data values.
#' @param x vector of arguments of prediction function.
#' @param type type of prediction function ("df" for density, "qf" for quantile function
#' or "dens" for density).
#'
#' @export
#'
setMethod("predict", c(object = "vtscopula"), function(object, data, x, type = "df") {
  if ((object@Vtransform@name != "Vlinear") | (!is(object@Wcopula,"swncopula")))
    stop("Fast forecasting restricted to linear V-transform and SWN Wcopula")
  delta <- object@Vtransform@pars[1]
  if (sum(x == delta) > 0)
    stop("points at fulcrum: change values")
  Udata <- as.numeric(vtrans(object@Vtransform, data))
  switch(type,
         "df" = {
           mult <- ifelse(x <= delta, -delta, 1-delta)
           ux <- vtrans(object@Vtransform, x)
           delta + mult * predict(object@Vcopula, Udata, ux, "df")
         },
         "qf" = {
           arg <- ifelse(x <= delta, (delta-x)/delta, (x-delta)/(1-delta))
           v <- predict(object@Vcopula, Udata, arg, "qf")
           u <- vinverse(object@Vtransform, v)
           ifelse(x <= delta, u, u+v)
         },
         "dens" = {
           ux <- vtrans(object@Vtransform, x)
           predict(object@Vcopula, Udata, ux, "dens")
         })
})


#' Extract parameters of vtscopula
#'
#' @param object an object of class \linkS4class{vtscopula}.
#'
#' @return A list of parameters.
#' @export
#' @keywords internal
#'
vtparlist <- function(object) {
    output <- object@Vcopula@pars
    vpars <- coef(object@Vtransform)
    if (length(vpars) > 1){
      vpars <- vpars[-1]
      output$vt <- vpars
    }
    if (!is(object@Wcopula, "swncopula")) {
      output$wcopula <- object@Wcopula@pars[[1]]
    }
    output
}

#' @describeIn vtscopula Simulation method for vtscopula class
#'
#' @param object an object of the class.
#' @param n length of realization.
#'
#' @export
#'
#' @examples
#' copobject <- armacopula(pars = list(ar = 0.6, ma = 0.2))
#' sim(vtscopula(copobject, Vtransform = V2p()))
setMethod("sim", c(object = "vtscopula"), function(object, n = 1000) {
  U <- sim(object@Vcopula, n)
  stochinverse(object@Vtransform, U, object@Wcopula)
})

#' Fit method for vtscopula class
#'
#' Fit object of class \linkS4class{vtscopula}
#' to data using maximum likelihood.
#'
#' @param x an object of class \linkS4class{vtscopula}.
#' @param y a vector or time series of data.
#' @param tsoptions list of optional arguments:
#' hessian is logical variable specifying whether Hessian matrix should be returned;
#' method is choice of optimization method.
#' @param control list of control parameters to be passed to the
#' \code{\link[stats]{optim}} function.
#'
#'
#' @return An object of class \linkS4class{tscopulafit}.
#' @export
#'
#' @examples
#' copobject <- armacopula(pars = list(ar = 0.6, ma = 0.2))
#' vtcop <- vtscopula(copobject, Vtransform = V2p())
#' y <- sim(vtcop)
#' fit(vtcop, y)
setMethod(
  "fit", c(x = "vtscopula", y = "ANY"),
  function(x, y,
           tsoptions = list(),
           control = list(maxit = 2000, warn.1d.NelderMead = FALSE)) {
    defaults <- list(hessian = FALSE, method = "Nelder-Mead")
    tsoptions <- setoptions(tsoptions, defaults)
    if (is(x, "tscopulafit")) {
      x <- x@tscopula
    }
    fulcrum <- as.numeric(x@Vtransform@pars["delta"])
    if (is.na(fulcrum))
      stop("V-transform must contain a fulcrum value")
    if (sum(y == fulcrum) > 0)
      stop("Data point at fulcrum")
    parlist <- vtparlist(x)
    fit <- optim(
      par = unlist(parlist),
      fulcrum = fulcrum,
      fn = vtscopula_objective,
      modelspec = x@Vcopula@modelspec,
      modeltype = is(x@Vcopula)[[1]],
      vt = x@Vtransform,
      wcopula = setwcopula(x),
      y = as.numeric(y),
      method = tsoptions$method,
      hessian = tsoptions$hessian,
      control = control
    )
    newpars <- relist(fit$par, parlist)
    newpars$vt <- c(delta = fulcrum, newpars$vt)
    x@Vtransform@pars <- newpars$vt
    if (!is(x@Wcopula, "swncopula")) {
      x@Wcopula@pars[[1]] <- newpars$wcopula
    }
    x@Vcopula@pars <- newpars[(names(newpars) != "vt") & (names(newpars) != "wcopula")]
    new("tscopulafit", tscopula = x, data = y, fit = fit)
  }
)

#' Extract W-copula
#'
#' @param x an object of class \linkS4class{tscopula}.
#'
#' @return A description of the W-copula.
#' @export
#' @keywords internal
#'
setwcopula <- function(x) {
  if (is(x@Wcopula, "swncopula")) {
    return(NULL)
  } else if (length(x@Wcopula@modelspec) > 1) {
    stop("W copula must be a pair copula: recommended is Frank")
  } else {
    return(x@Wcopula@modelspec[[1]])
  }
}

#' Objective function for vtscopula fitting
#'
#' @param theta vector of parameters
#' @param fulcrum fixed value for fulcrum
#' @param modelspec list containing tscopula specification
#' @param modeltype character vector specifying tscopula type
#' @param vt a v transform
#' @param wcopula W copula for serial dependence
#' @param y vector of data
#'
#' @return Value of objective function at parameters.
#' @export
#' @keywords internal
#'
vtscopula_objective <- function(theta, fulcrum, modelspec, modeltype, vt, wcopula, y) {
  n_corepars <- switch(modeltype,
    armacopula = sum(modelspec),
    dvinecopula = sum(sapply(modelspec, function(v) {
      v$npars
    })),
    dvinecopulaNE = sum(sapply(modelspec, function(v) {
      v$npars
    })),
    dvinecopula2 = modelspec$npar,
    dvinecopula3 = modelspec$npar
  )
  theta_core <- theta[1:n_corepars]
  theta_vt <- c(delta = fulcrum)
  n_vtextra <- length(theta[substring(names(theta), 1, 2) == "vt"])
  if (n_vtextra > 0) {
    theta_vtextra <- theta[substring(names(theta), 1, 2) == "vt"]
    names(theta_vtextra) <- substring(names(theta_vtextra), 4)
    theta_vt <- c(theta_vt, theta_vtextra)
    }
  if ((min(theta_vt) < 0) | (theta_vt["delta"] > 1))
    return(NA)
  V <- do.call(vt@Vtrans, append(theta_vt, list(u = y)))
  if (min(V) == 0)
    return(NA)
  objective <- eval(parse(text = paste(modeltype, "_objective", sep = "")))
  value <- objective(theta_core, modelspec, V)
  if (!(is.null(wcopula))) {
    theta_wcopula <- theta[substring(names(theta), 1, 7) == "wcopula"]
    value <- value + wobjective(theta_wcopula, wcopula, theta_vt, vt@gradient, y)
  }
  value
}

#' Additional objective for generalized processes
#'
#' @param theta parameter of W copula
#' @param paircop specification of W copula
#' @param theta_vt parameter(s) of v-transform
#' @param vtgrad gradient of v-transform
#' @param u vector of data
#'
#' @return Value of additional objective for W-copula.
#' @keywords internal
#'
wobjective <- function(theta, paircop, theta_vt, vtgrad, u) {
  if (length(theta_vt) > 0) {
    delta <- theta_vt["delta"]
  } else {
    delta <- 0.5
  }
  if (length(theta_vt) <= 1) {
    Delta <- rep(delta, length(u))
  } else {
    Delta <- -1 / do.call(vtgrad, append(theta_vt, list(u = u)))
    Delta[u > delta] <- Delta[u > delta] + 1
  }
  Delta1 <- Delta[-length(Delta)]
  Delta2 <- Delta[-1]
  u1 <- u[-length(u)]
  u2 <- u[-1]
  zone1 <- (u1 <= delta) & (u2 <= delta)
  zone2 <- (u1 <= delta) & (u2 > delta)
  zone3 <- (u1 > delta) & (u2 <= delta)
  zone4 <- (u1 > delta) & (u2 > delta)
  pars <- theta[1]
  if (paircop$npars > 1) {
    pars <- c(pars, theta[2])
  }
  copfamily <- tryCatch(rvinecopulib::bicop_dist(
    family = tolower(paircop$family),
    rotation = paircop$rotation,
    parameters = pars
  ),
  error = function(e) {
    return(NA)
  }
  )
  if (is.na(copfamily[[1]])) {
    return(NA)
  }
  Cvals <- rvinecopulib::pbicop(u = cbind(Delta1, Delta2), family = copfamily)
  output <- rep(NA, length(u) - 1)
  output[zone1] <- Cvals[zone1] / (Delta1[zone1] * Delta2[zone1])
  output[zone2] <- (Delta1[zone2] - Cvals[zone2]) / (Delta1[zone2] * (1 - Delta2[zone2]))
  output[zone3] <- (Delta2[zone3] - Cvals[zone3]) / ((1 - Delta1[zone3]) * Delta2[zone3])
  output[zone4] <- (1 - Delta1[zone4] - Delta2[zone4] + Cvals[zone4]) / ((1 - Delta1[zone4]) * (1 - Delta2[zone4]))
  -sum(log(output))
}

#' Profile likelihood for fulcrum parameter
#'
#' @param data a vector or time series of data on (0,1).
#' @param tscopula an object of class \linkS4class{tscopulaU} or \linkS4class{vtscopula}.
#' @param locations vector containing locations of different values for fulcrum.
#' @param plot logical values specifying whether plot should be created.
#'
#' @return A matrix containing fulcrum values and log likelihood values.
#' @export
#'
#' @examples
#' copobject <- armacopula(pars = list(ar = 0.6, ma = 0.2))
#' vtcop <- vtscopula(copobject, Vtransform = V2p())
#' y <- sim(vtcop)
#' profilefulcrum(y, vtcop)
profilefulcrum <- function(data,
                           tscopula = dvinecopula(family = 1, pars = list(0.1)),
                           locations = seq(0, 1, by = 0.1),
                           plot = TRUE) {
  if (!is(tscopula, "vtscopula")) {
    tscopula <- new("vtscopula",
      Vcopula = tscopula,
      Vtransform = Vlinear(),
      Wcopula = swncopula()
    )
  }
  if (length(tscopula@Vtransform@pars) == 0) {
    tscopula@Vtransform <- Vlinear()
  }
  results <- numeric(length(locations))
  for (i in seq_along(locations)) {
      coptofit <- tscopula
    if (locations[i] %in% c(0,1))
      coptofit@Vtransform <- Vlinear()
    coptofit@Vtransform@pars["delta"] <- locations[i]
    fitted_model <- fit(
      x = coptofit,
      y = data,
      tsoptions = list(
        hessian = FALSE
      )
    )
    results[i] <- -fitted_model@fit$value
  }
  if (plot) {
    plot(locations, results, xlab = expression(delta), ylab = "L", type = "l")
  }
  invisible(cbind(fulcrum = locations, logLik = results))
}

#' Fit tscm Jointly
#'
#' @param x an object of class \linkS4class{tscm}.
#' @param y a vector or time series of data.
#' @param tsoptions list of variables
#' @param control list of control parameters to be passed to the
#' \code{\link[stats]{optim}} function.
#'
#' @return An object of class \linkS4class{tscmfit}.
#' @keywords internal
#'
fitFULLb <- function(x, y, tsoptions, control) {
  dens <- eval(parse(text = paste("d", x@margin@name, sep = "")))
  cdf <- eval(parse(text = paste("p", x@margin@name, sep = "")))
  parlist <- vtparlist(x@tscopula)
  parlist$margin <- x@margin@pars
  fulcrum <- as.numeric(x@tscopula@Vtransform@pars["delta"])
  if (tsoptions$changeatzero){
    if (length(y[y == 0]) > 0)
      stop("Remove zeros in dataset")
    fulcrum <- NA
  }
  fit <- optim(
    par = unlist(parlist),
    fulcrum = fulcrum,
    fn = tsc_objectiveb,
    modelspec = x@tscopula@Vcopula@modelspec,
    modeltype = is(x@tscopula@Vcopula)[[1]],
    dens = dens,
    cdf = cdf,
    vt = x@tscopula@Vtransform,
    wcopula = setwcopula(x@tscopula),
    y = as.numeric(y),
    method = tsoptions$method,
    hessian = tsoptions$hessian,
    control = control
  )
  newpars <- relist(fit$par, parlist)
  x@margin@pars <- newpars$margin
  newpars <- newpars[names(newpars) != "margin"]
  if (is.na(fulcrum))
    fulcrum <- pmarg(x@margin, 0)
  newpars$vt <- c(delta = fulcrum, newpars$vt)
  x@tscopula@Vtransform@pars <- newpars$vt
  newpars <- newpars[names(newpars) != "vt"]
  if (!is(x@tscopula@Wcopula, "swncopula")) {
    x@tscopula@Wcopula@pars[[1]] <- newpars$wcopula
    newpars <- newpars[names(newpars) != "wcopula"]
  }
  x@tscopula@Vcopula@pars <- newpars
  new("tscmfit", tscopula = x@tscopula, margin = x@margin, data = y, fit = fit)
}

#' Objective function for full fit with v-transform
#'
#' @param theta vector of parameters
#' @param fulcrum value for fulcrum
#' @param modelspec list containing model specification
#' @param modeltype character vector giving type of model
#' @param dens marginal density function
#' @param cdf marginal distribution function
#' @param vt v-transform
#' @param wcopula W copula if used
#' @param y vector of data
#'
#' @return Value of objective function at parameters.
#' @keywords internal
#'
tsc_objectiveb <-
  function(theta, fulcrum, modelspec, modeltype, dens, cdf, vt, wcopula, y) {
    margpars <- theta[substring(names(theta), 1, 6) == "margin"]
    nonmargpars <- theta[substring(names(theta), 1, 6) != "margin"]
    names(margpars) <- substring(names(margpars), 8)
    dx <- do.call(dens, append(as.list(margpars), list(x = y, log = TRUE)))
    termA <- -sum(dx)
    if (is.na(termA)) {
      return(NA)
    }
    U <- do.call(cdf, append(margpars, list(q = y)))
    if (is.na(fulcrum))
      fulcrum <- do.call(cdf, append(margpars, list(q = 0)))
    termBC <-
      vtscopula_objective(nonmargpars, fulcrum, modelspec, modeltype, vt, wcopula, U)
    return(termA + termBC)
  }

#' @describeIn vtscopula Calculate Kendall's tau values for vtscopula model
#'
#' @param object an object of the class.
#' @param lagmax maximum value of lag.
#'
#' @export
#'
#' @examples
#' mod <- vtscopula(armacopula(list(ar = 0.95, ma = -0.85)))
#' kendall(mod)
setMethod("kendall", c(object = "vtscopula"), function(object, lagmax = 20){
kendall(object@Vcopula, lagmax)
}
)
