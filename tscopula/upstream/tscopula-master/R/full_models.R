#' Full models
#'
#' Class of objects for composite time series models consisting
#' of stationary copula processes and marginal distributions.
#'
#' @slot tscopula an object of class \linkS4class{tscopula}.
#' @slot margin an object of class \linkS4class{margin}.
#'
#' @export
#'
setClass("tscm",
  slots = list(
    tscopula = "tscopula",
    margin = "margin"
  )
)

#' @describeIn tscm Show method for tscm class
#'
#' @param object an object of the class.
#'
#' @export
#'
setMethod("show", "tscm", function(object) {
  cat("object class: ", is(object)[[1]], "\n", sep = "")
  cat("_______ \n")
  cat("MARGIN: \n")
  show(object@margin)
  cat("_______\n")
  cat("COPULA: \n")
  show(object@tscopula)
  if (is(object, "tscmfit")) {
    ests <- object@fit$par
    if ("skeleton" %in% names(attributes(ests))) {
      attributes(ests)$skeleton <- NULL
    }
    if (is.element("hessian", names(object@fit))) {
      ses <- safe_ses(object@fit$hessian)
      ests <- rbind(ests, ses)
      dimnames(ests)[[1]] <- c("par", "se")
    }
    cat("_________________________\n")
    cat("summary of all estimates:\n")
    print(ests)
    cat("convergence status:", object@fit$convergence, ", log-likelihood:",
      -object@fit$value, "\n",
      sep = " "
    )
  }
})

#' @describeIn tscm Coefficient method for tscm class
#'
#' @param object an object of the class.
#'
#' @export
#'
setMethod("coef", "tscm", function(object) {
  c(coef(object@tscopula), coef(object@margin))
})

#' Constructor function for time series
#'
#' @param tscopula an object of class \linkS4class{tscopula}.
#' @param margin an object of class \linkS4class{margin}.
#'
#' @return An object of class \linkS4class{tscm}.
#' @export
#'
#' @examples
#' tscm(dvinecopula(family = "gauss", pars = 0.5), margin("doubleweibull"))
tscm <- function(tscopula, margin = new("margin", name = "unif")) {
  if (is(tscopula, "tscopulafit")) {
    tscopula <- tscopula@tscopula
  }
  if (is(margin, "marginfit")) {
    margin <- margin@margin
  }
  new("tscm",
    tscopula = tscopula,
    margin = margin
  )
}

#' @describeIn tscm Simulation method for tscm class
#'
#' @param object an object of the class.
#' @param n length of realization.
#'
#' @export
#'
#' @examples
#' mod <- tscm(dvinecopula(family = "gauss", pars = 0.5), margin("doubleweibull"))
#' sim(mod)
setMethod(
  "sim",
  c(object = "tscm"),
  function(object, n = 1000) {
    Utilde <- sim(object@tscopula, n)
    qmarg(object@margin, Utilde)
  }
)

#' @describeIn tscm Prediction method for tscm class
#'
#' @param object an object of the class.
#' @param data vector of past data values.
#' @param x vector of arguments of prediction function.
#' @param type type of prediction function ("df" for density, "qf" for quantile function
#' or "dens" for density).
#' @param qtype type of empirical quantile estimate.
#' @param proper logical variable stating whether the standard empirical distribution function should
#' be used when the margin is empirical; otherwise an improper distribution that is bounded away from 0
#' and 1 is used.
#'
#' @export
#'
setMethod("predict", c(object = "tscm"), function(object, data, x, type = "df", qtype = 7, proper = FALSE) {
  margmod <- object@margin
  if(margmod@name == "edf")
    return(predict_empirical(object, data, x, type, qtype, proper))
  Udata <- pmarg(margmod, data)
  uvals <- switch(type,
                  "df" = pmarg(margmod, x),
                  "qf" = x,
                  "dens" = pmarg(margmod, x))
  upred <- predict(object@tscopula, Udata, uvals, type = type)
  switch(type,
         "df" = upred,
         "qf" = qmarg(margmod, upred),
         "dens" = upred * dmarg(margmod, x))
})

#' @describeIn tscm Calculate Kendall's tau values for pair copulas for tscm class
#'
#' @param object an object of the class.
#' @param lagmax maximum value of lag.
#'
#' @export
#'
setMethod("kendall", c(object = "tscm"), function(object, lagmax = 20) {
  kendall(object@tscopula, lagmax)
}
)

#' Adjusted empirical distribution function
#'
#' @param x argument of empirical distribution function.
#' @param data vector of data for constructing empirical
#' distribution function.
#' @param proper logical variable which when set to TRUE will return
#' the standard empirical distribution function.
#'
#' @return a vector of same length as x
#' @export
#'
pedf <- function(x, data, proper = FALSE){
  df <- ecdf(data)(x)
  if (!proper){
    n <- length(data)
    df <- df*n/(n+1)
    df[x < min(data)] <- 0.5/(n+1)
    df[x > max(data)] <- (n+0.5)/(n+1)
  }
  df
}

#' Prediction function for tscm class with empirical margin
#'
#' @param object an object of the class.
#' @param data vector of past data values.
#' @param x vector of arguments of prediction function.
#' @param type type of prediction function ("df" for density, "qf" for quantile function
#' or "dens" for density).
#' @param qtype type of empirical quantile estimate.
#' @param proper logical variable stating whether the standard empirical distribution function should
#' be used when the margin is empirical; otherwise an improper distribution that is bounded away from 0
#' and 1 is used.
#' @keywords internal
#'
predict_empirical <- function(object, data, x, type, qtype, proper){
  Udata <- strank(data)
  uvals <- switch(type,
                  "df" = pedf(x, data, proper),
                  "qf" = x,
                  "dens" = pedf(x, data, proper))
  upred <- predict(object@tscopula, Udata, uvals, type = type)
  switch(type,
         "df" = upred,
         "qf" = quantile(data, upred, type = qtype),
         "dens" = upred * kdensity::kdensity(data)(x))
}

#' Fitted tscm model
#'
#' Class of objects for fitted \linkS4class{tscm} models.
#'
#' @slot tscopula an object of class \linkS4class{tscopula}.
#' @slot margin an object of class \linkS4class{margin}.
#' @slot data a vector or time series of data to which process has been fitted.
#' @slot fit a list containing details of the fit.
#'
#' @export
#'
setClass(
  "tscmfit",
  contains = "tscm",
  slots = list(
    tscopula = "tscopula",
    margin = "margin",
    data = "ANY",
    fit = "list"
  )
)

#' Fit method for tscm class
#'
#' @param x an object of class \linkS4class{tscm}.
#' @param y a vector or time series of data.
#' @param tsoptions a list of parameters passed to fitting.
#' @param control list of control parameters to be passed to the
#' \code{\link[stats]{optim}} function.
#' @param method character string specifying method.
#'
#'
#' @return An object of class \linkS4class{tscmfit}.
#' @export
#'
#' @examples
#' mod <- tscm(dvinecopula(family = "gauss", pars = 0.5), margin("doubleweibull"))
#' y <- sim(mod)
#' fit(mod, y)
setMethod(
  "fit", c(x = "tscm", y = "ANY"),
  function(x,
           y,
           tsoptions = list(),
           control = list(),
           method = "IFM") {
    defaults <- list(hessian = FALSE, method = "Nelder-Mead", changeatzero = FALSE)
    cdefaults <- list(warn.1d.NelderMead = FALSE, trace = FALSE, maxit = 5000)
    tsoptions <- setoptions(tsoptions, defaults)
    control <- setoptions(control, cdefaults)
    if (is(x, "tscmfit")) {
      tscopula <- x@tscopula
      margin <- x@margin
      if (is(tscopula, "tscopulafit")) {
        tscopula <- tscopula@tscopula
      }
      if (is(margin, "marginfit")) {
        margin <- margin@margin
      }
      x <- new("tscm", tscopula = tscopula, margin = margin)
    }
    if (x@margin@name == "edf")
      method <- "empirical"
    if ((method == "full") & is(x@tscopula, "tscopulaU")) {
      method <- paste(method, "A", sep = "")
    }
    if ((method == "full") & is(x@tscopula, "vtscopula")) {
      method <- paste(method, "B", sep = "")
    }
    switch(method,
      empirical = fitEDF(x, y, tsoptions, control),
      IFM = fitSTEPS(x, y, tsoptions, control),
      fullA = fitFULLa(x, y, tsoptions, control),
      fullB = fitFULLb(x, y, tsoptions, control),
      stop("Not a known method")
    )
  }
)

#' Fit tscm using empirical distribution function
#'
#' @param x an object of class \linkS4class{tscm}.
#' @param y a vector or time series of data.
#' @param tsoptions a list of parameters passed to fitting.
#' This differs according to the class of x.
#' @param control list of control parameters to be passed to the
#' \code{\link[stats]{optim}} function.
#'
#' @return An object of class \linkS4class{tscmfit}.
#' @keywords internal
#'
fitEDF <- function(x, y, tsoptions, control) {
  U <- strank(as.numeric(y))
  if (is(x@tscopula, "vtscopula") & tsoptions$changeatzero){
    if (length(y[y == 0]) > 0)
      stop("Remove zeros in dataset")
    x@tscopula@Vtransform@pars["delta"] <-
      as.numeric((length(y[y < 0]) + 0.5)/(length(y) + 1))
  }
  copfit <- fit(x@tscopula, U, tsoptions, control)
  new("tscmfit",
    tscopula = copfit@tscopula,
    margin = new("margin", name = "edf"),
    data = y,
    fit = copfit@fit
  )
}

#' Fit tscm in two steps
#'
#' @param x an object of class \linkS4class{tscm}.
#' @param y a vector or time series of data.
#' @param tsoptions a list of parameters passed to fitting.
#' This differs according to the class of x.
#' @param control list of control parameters to be passed to the
#' \code{\link[stats]{optim}} function.
#'
#' @return An object of class \linkS4class{tscmfit}.
#' @keywords internal
#'
fitSTEPS <- function(x, y, tsoptions, control) {
  margfit <- fit(x@margin, y, tsoptions = tsoptions, control = control)
  U <- pmarg(margfit, y)
  if (is(x@tscopula, "vtscopula") & tsoptions$changeatzero){
    if (length(y[y == 0]) > 0)
      stop("Remove zeros in dataset")
    x@tscopula@Vtransform@pars["delta"] <- as.numeric(pmarg(margfit, 0))
  }
  copfit <- fit(x@tscopula, U, tsoptions = tsoptions, control = control)
  combinedfit <- list()
  names(margfit@fit$par) <- paste("margin.", names(margfit@fit$par), sep = "")
  combinedfit$par <- c(copfit@fit$par, margfit@fit$par)
  if ("hessian" %in% names(tsoptions)) {
    if (tsoptions$hessian) {
      combinedfit$hessian <- Matrix::bdiag(margfit@fit$hessian, copfit@fit$hessian)
    }
  }
  combinedfit$convergence <- sum(copfit@fit$convergence, margfit@fit$convergence)
  combinedfit$value <- sum(copfit@fit$value, margfit@fit$value)
  new("tscmfit",
    tscopula = copfit@tscopula,
    margin = margfit@margin,
    data = y,
    fit = combinedfit
  )
}

#' Fit tscm jointly
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
fitFULLa <- function(x, y, tsoptions, control) {
  dens <- eval(parse(text = paste("d", x@margin@name, sep = "")))
  cdf <- eval(parse(text = paste("p", x@margin@name, sep = "")))
  parlist <- x@tscopula@pars
  parlist$margin <- x@margin@pars
  fit <- optim(
    par = unlist(parlist),
    fn = tsc_objectivea,
    modelspec = x@tscopula@modelspec,
    modeltype = is(x@tscopula)[[1]],
    dens = dens,
    cdf = cdf,
    y = as.numeric(y),
    method = tsoptions$method,
    hessian = tsoptions$hessian,
    control = control
  )
  newpars <- relist(fit$par, parlist)
  x@margin@pars <- newpars$margin
  x@tscopula@pars <- newpars[names(newpars) != "margin"]
  new("tscmfit", tscopula = x@tscopula, margin = x@margin, data = y, fit = fit)
}



#' Objective function for full of tscopula plus margin model
#'
#' @param theta vector of parameter values
#' @param modelspec list containing model specification
#' @param modeltype character string giving type of model
#' @param dens marginal density function
#' @param cdf marginal cdf
#' @param y vector of data values
#'
#' @return Value of objective function at parameters.
#' @keywords internal
#'
tsc_objectivea <-
  function(theta, modelspec, modeltype, dens, cdf, y) {
    margpars <- theta[substring(names(theta), 1, 6) == "margin"]
    nonmargpars <- theta[substring(names(theta), 1, 6) != "margin"]
    names(margpars) <- substring(names(margpars), 8)
    dx <- do.call(dens, append(as.list(margpars), list(x = y, log = TRUE)))
    termA <- -sum(dx)
    if (is.na(termA)) {
      return(NA)
    }
    U <- do.call(cdf, append(margpars, list(q = y)))
    objective <- eval(parse(text = paste(modeltype, "_objective", sep = "")))
    termBC <-
      objective(nonmargpars, modelspec, U)
    return(termA + termBC)
  }

#' Convert tscopula object to tscm object
#'
#' @param from a \linkS4class{tscopula} object.
#' @param to a \linkS4class{tscm} object.
#' @param strict logical variable stating whether strict coercion should be enforced.
#'
#' @return A \linkS4class{tscm} object.
#' @export
#'
setMethod(
  "coerce", c(from = "tscopula", to = "tscm"),
  function(from, to = "tsc", strict = TRUE) {
    new("tscm", tscopula = from, margin = new("margin", name = "unif"))
  }
)

#' Convert tscopulafit object to be tscmfit object
#'
#' @param from a \linkS4class{tscopulafit} object.
#' @param to a \linkS4class{tscmfit} object.
#' @param strict logical variable stating whether strict coercion should be enforced.
#'
#' @return A \linkS4class{tscmfit} object.
#' @export
#'
setMethod(
  "coerce", c(from = "tscopulafit", to = "tscmfit"),
  function(from, to = "tscmfit", strict = TRUE) {
    new("tscmfit",
      tscopula = from,
      margin = new("margin", name = "unif"),
      data = from@data,
      fit = from@fit
    )
  }
)

#' @describeIn tscmfit method for tscmfit class
#'
#' @param object an object of the class.
#'
#' @export
#'
setMethod("logLik", "tscmfit", function(object) {
  ll <- -object@fit$value[length(object@fit$value)]
  attr(ll, "nobs") <- length(object@data)
  attr(ll, "df") <- length(object@fit$par)
  class(ll) <- "logLik"
  ll
})

#' @describeIn tscmfit Residual method for tscmfit class
#'
#' @param object an object of the class.
#' @param trace extract trace instead of residuals.
#'
#' @export
#'
setMethod("resid", "tscmfit",
          function(object, trace = FALSE) {
            object <- new("tscopulafit",
                tscopula = object@tscopula,
                data = pmarg(object@margin, object@data),
                fit = object@fit)
            resid(object, trace)
          })

#' @describeIn tscmfit Prediction method for tscmfit class
#'
#' @param object an object of the class.
#' @param x vector of arguments of prediction function.
#' @param type type of prediction function ("df" for density, "qf" for quantile function
#' or "dens" for density).
#' @param qtype type of empirical quantile estimate.
#' @param proper logical variable stating whether the standard empirical distribution function should
#' be used when the margin is empirical; otherwise an improper distribution that is bounded away from 0
#' and 1 is used.
#'
#' @export
#'
setMethod("predict", c(object = "tscmfit"), function(object, x, type = "df", qtype = 7, proper = FALSE) {
  predict(tscm(object@tscopula, object@margin), object@data, x, type = type, qtype = qtype, proper = proper)
})

#' Plot method for tscmfit class
#'
#' @param x an object of class \linkS4class{tscmfit}.
#' @param plottype type of plot required.
#' @param bw logical variable specifying whether black-white options should be chosen.
#' @param lagmax maximum lag value for dvinecopula2 plots
#'
#' @return No return value, generates plot.
#' @export
#'
setMethod("plot", c(x = "tscmfit", y = "missing"),
          function(x, plottype = "residual", bw = FALSE, lagmax = 30) {
            if (x@margin@name == "edf")
              stop("These plots require parametric margins")
            if (plottype == "margin"){
              marginmod <- new("marginfit", margin = x@margin, data = x@data, fit = x@fit)
              plot(marginmod, bw = bw)
            }
            else if (plottype == "volproxy")
              plot_volproxy(x, bw = bw)
            else if (plottype == "volprofile")
              plot_volprofile(x, bw = bw)
            else if (plottype %in% c("residual", "kendall", "glag", "vtransform")){
              copmod <- new("tscopulafit", tscopula = x@tscopula, data = pmarg(x@margin, x@data), fit = list(NULL))
              plot(copmod, plottype = plottype, bw = bw, lagmax = lagmax)
            }
            else
              stop("Not a valid plot type")
          })


#' Plot function for volatility profile plot
#'
#' @param x an object of class \linkS4class{tscmfit}.
#' @param bw logical variable specifying whether black-white options should be chosen.
#'
#' @return No return value, generates plot.
#' @keywords internal
#'
plot_volprofile <- function(x, bw) {
  if (!(is(x@tscopula, "vtscopula")))
    stop("tscopula must be vtscopula")
  xvals <- seq(from = 0,
               to = max(abs(x@data)),
               length = 100)
  vpars <- x@tscopula@Vtransform@pars
  brk <- qmarg(x@margin, vpars["delta"])
  u <- pmarg(x@margin, brk - xvals)
  v <- vtrans(x@tscopula@Vtransform, u)
  yvals <- qmarg(x@margin, u + v) - brk
  plot(xvals, yvals, type = "l", ylab = expression(g[T](x)))
  colchoice <- ifelse(bw, "gray50", "red")
  abline(0, 1, col = colchoice, lty = 2)
}

#' Plot function for volatility proxy plot
#'
#' @param x an object of class \linkS4class{tscmfit}.
#' @param bw logical variable specifying whether black-white options should be chosen.
#'
#' @return No return value, generates plot.
#' @keywords internal
#'
plot_volproxy <- function(x, bw){
  if (!(is(x@tscopula, "vtscopula")))
    stop("tscopula must be vtscopula")
  X <- as.numeric(x@data)
  U <- pmarg(x@margin, X)
  V <- vtrans(x@tscopula@Vtransform, U)
  colchoice <- ifelse(bw, "gray50", "red")
  plot(X, qnorm(strank(V)), xlab = "data", ylab = "std. vol proxy", col = colchoice)
  lines(sort(X), qnorm(V[order(X)]))
}
