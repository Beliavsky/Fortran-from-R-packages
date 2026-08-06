#' @title Time series copulas of class tscopulaU
#'
#' @description S4 Class union for basic time series copula types. These are
#' \linkS4class{armacopula}, \linkS4class{dvinecopula} and \linkS4class{dvinecopula2},
#'
#' @exportClass tscopulaU
#'
setClassUnion("tscopulaU", c("armacopula", "sarmacopula",
                             "dvinecopula", "dvinecopula2", "dvinecopula3"))

#' Calculate standardized ranks of data
#'
#' @param x a vector or time series of data.
#'
#' @return A vector or time series of standardized ranks in the interval (0,1)
#' @export
#'
#' @examples
#' strank(rnorm(100))
strank <- function(x) {
  U <- rank(x, na.last = "keep", ties.method = "random") / (length(x) + 1)
  if (inherits(x, "zoo")) {
    attributes(U) <- attributes(x)
  }
  U
}

#' Generic for estimating time series models
#'
#' Methods are available for objects of class \linkS4class{tscopulaU},
#' \linkS4class{vtscopula}, \linkS4class{tscopulafit}, \linkS4class{margin} and
#' \linkS4class{tscm}.
#'
#' @param x an object of the model class.
#' @param y a vector or time series of data.
#' @param ... further arguments to be passed on.
#'
#' @return An object of the fitted model class.
#' @export
#'
#'
setGeneric("fit", function(x, y, ...) {
  standardGeneric("fit")
})

#' Fitted time series copula processes
#'
#' Class of objects for fitted time series copula processes.
#'
#' @slot tscopula an object of class \linkS4class{tscopula}.
#' @slot data a vector or time series of data.
#' @slot fit a list containing details of the fit.
#'
#' @export
#'
setClass("tscopulafit",
  contains = c("tscopula"),
  slots = list(
    tscopula = "tscopula",
    data = "ANY",
    fit = "list"
  )
)

#' @describeIn tscopulafit Simulation method for tscopulafit class
#'
#' @param object an object of class \linkS4class{tscopulafit}.
#' @param n length of realization.
#'
#' @export
#'
#' @examples
#' ar1 <- armacopula(list(ar = 0.7))
#' data <- sim(ar1, 1000)
#' ar1fit <- fit(ar1, data)
#' sim(ar1fit)
setMethod("sim", c(object = "tscopulafit"), function(object, n = 1000) {
  sim(object@tscopula, n)
})

#' @describeIn tscopulafit Calculate Kendall's tau values for pair copulas for tscopulafit class
#'
#' @param object an object of the class.
#' @param lagmax maximum value of lag.
#'
#' @export
#'
setMethod("kendall", c(object = "tscopulafit"), function(object, lagmax = 20) {
  kendall(object@tscopula, lagmax)
}
)


#' @describeIn tscopulafit Coef method for tscopulafit class
#'
#' @param object an object of class \linkS4class{tscopulafit}.
#'
#' @export
#'
setMethod("coef", c(object = "tscopulafit"), function(object) {
  coef(object@tscopula)
})

#' @describeIn tscopulafit Show method for tscopulafit objects
#'
#' @param object an object of class \linkS4class{tscopulafit}.
#'
#' @export
#'
setMethod("show", "tscopulafit", function(object) {
  if (is(object@tscopula, "vtscopula")) {
    show(object@tscopula)
  } else if (is(object@tscopula, "tscopulaU")) {
    cat("object class: ", is(object@tscopula)[[1]], "\n", sep = "")
    cat("name: ", object@tscopula@name, "\n", sep = "")
    if (is(object@tscopula, "dvinecopula")) {
      fams <- sapply(object@tscopula@modelspec, FUN = function(v) {
        tmp <- v$family
        if (v$rotation != 0) {
          tmp <- paste(tmp, v$rotation, sep = "")
        }
        if (is(object, "dvinecopulaNE")) {
          if (!(v$exchangeable)) {
            tmp <- paste(tmp, "_NE", sep = "")
          }
        }
        tmp
      })
      if (length(unique(fams)) == 1) {
        fams <- fams[1]
      }
      cat("copula family: ", fams, "\n")
    }
    if (is(object@tscopula, "dvinecopula2")) {
      modobject <- object@tscopula
      famname <- modobject@modelspec$family
      if (modobject@modelspec$rotation !=0)
        famname <- paste(famname, "with rotation", modobject@modelspec$rotation)
      cat("copula family: ", famname, "\n", sep = "")
      if (modobject@modelspec$negtau != "none")
        cat("negative tau treatment: ", modobject@modelspec$negtau, "\n", sep = "")
      kpacf  <- modobject@modelspec$kpacf
      if (modobject@modelspec$maxlag != Inf)
        kpacf <- paste(kpacf, "with max lag", modobject@modelspec$maxlag)
      cat("KPACF: ", kpacf,"\n", sep = "")
    }
    if (is(object@tscopula, "dvinecopula3")) {
      modobject <- object@tscopula
      cat("non-Gaussian substitutions:\n")
      locs <- modobject@modelspec$location
      cat(" - locations:", locs, "\n", sep = " ")
      cat(" - families:", modobject@modelspec$family, "\n", sep = " ")
      tau <- kendall(modobject, lagmax = max(locs))[locs]
      rot <- ifelse(tau >= 0, modobject@modelspec$posrot, modobject@modelspec$negrot)
      cat(" - rotations:", rot, "\n", sep = " ")
      cat(" - Kendall's tau:", round(tau,3), "\n", sep = " ")
      kpacf  <- modobject@modelspec$kpacf
      if (modobject@modelspec$maxlag != Inf)
        kpacf <- paste(kpacf, "with max lag", modobject@modelspec$maxlag)
      cat("KPACF: ", kpacf,"\n", sep = "")
    }
  }
  else {
    stop("Unknown tscopula type")
  }
  ests <- object@fit$par
  attributes(ests)$skeleton <- NULL
  if (is.element("hessian", names(object@fit))) {
    ses <- safe_ses(object@fit$hessian)
    ests <- rbind(ests, ses)
    dimnames(ests)[[1]] <- c("par", "se")
  }
  cat("_____________________\n")
  cat("Summary of estimates:\n")
  print(ests)
  cat("convergence status: ", object@fit$convergence, ", log-likelihood: ",
    -object@fit$value, "\n",
    sep = ""
  )
})

#' Calculate standard errors safely
#'
#' @param hess a Hessian matrix from a model fit.
#'
#' @return a vector of standard errors.
#' @export
#'
safe_ses <- function(hess) {
  hessinverse <- tryCatch(solve(hess), error = function(e) {
    warning("Hessian can't be inverted")
    return(diag(rep(NA, length(diag(hess)))))
  })
  sqrt(abs(diag(hessinverse)))
}

#' Fit method for tscopulafit class
#'
#' @param x an object of class \linkS4class{tscopulafit}.
#' @param y vector or time series of data to which the copula process is to be fitted.
#' @param tsoptions list of options
#' @param control list of control parameters to be passed to the
#' \code{\link[stats]{optim}} function.
#'
#' @return An object of class \linkS4class{tscopulafit}.
#' @export
#'
#' @examples
#' ar1 <- armacopula(list(ar = 0.7))
#' data <- sim(ar1, 1000)
#' ar1fit <- fit(fit(ar1, data), sim(ar1, 1000))
setMethod(
  "fit", c(x = "tscopulafit", y = "ANY"),
  function(x, y,
           tsoptions = list(),
           control = list(warn.1d.NelderMead = FALSE)) {
    x <- x@tscopula
    fit(x, y, tsoptions = tsoptions, control = control)
  }
)

#' Fit method for tscopulaU class
#'
#' @param x an object of class \linkS4class{tscopulaU}.
#' @param y vector or time series of data to which the copula process is to be fitted.
#' @param tsoptions list of options
#' @param control list of control parameters to be passed to the
#' \code{\link[stats]{optim}} function.
#'
#' @return An object of class \linkS4class{tscopulafit}.
#' @export
#'
#' @examples
#' data <- sim(armacopula(list(ar = 0.5, ma = 0.4)), n = 1000)
#' fit(armacopula(list(ar = 0.5, ma = 0.4)), data)
setMethod(
  "fit", c(x = "tscopulaU", y = "ANY"),
  function(x, y,
           tsoptions = list(),
           control = list()) {
    defaults <- list(hessian = FALSE, method = "Nelder-Mead")
    cdefaults <- list(warn.1d.NelderMead = FALSE, maxit = 1000)
    tsoptions <- setoptions(tsoptions, defaults)
    control <- setoptions(control, cdefaults)
    objective <- eval(parse(text = paste(is(x)[[1]], "_objective", sep = "")))
    fit <- optim(
      par = unlist(x@pars),
      fn = objective,
      modelspec = x@modelspec,
      u = as.numeric(y),
      method = tsoptions$method,
      hessian = tsoptions$hessian,
      control = control
    )
    x@pars <- relist(fit$par, x@pars)
    new("tscopulafit",
      tscopula = x,
      data = y,
      fit = fit
    )
  }
)

#' Set optional choices for tscopula fitting
#'
#' @param tsoptions a list of options chosen by user.
#' @param defaults a list of defaults specified in function.
#'
#' @return the definitive list of options to be used in function.
#' @keywords internal
#'
setoptions <- function(tsoptions, defaults) {
  defnames <- names(defaults)
  defaults[(optnames <- names(tsoptions))] <- tsoptions
  defaults
}

#' @describeIn tscopulafit logLik method for tscopulafit class
#'
#' @param object an object of class \linkS4class{tscopulafit}.
#'
#' @export
#'
setMethod("logLik", "tscopulafit", function(object) {
  logLik(as(object, "tscmfit"))
})

#' Akaike Corrected Information Criterion
#'
#' @param object a fitted model object for which there exists a logLik method
#' to extract the corresponding log-likelihood.
#' @param ... optionally more fitted model objects.
#'
#' @return If just one object is provided, a numeric value with the corresponding AICC value.
#'
#' If multiple objects are provided, a data.frame with rows corresponding to the objects and
#' columns representing the number of parameters in the model (df) and the AICC.
#' @export
#'
AICc <- function(object, ...)
{
  ll   <- if(isNamespaceLoaded("stats4")) stats4::logLik else logLik
  Nobs <- if(isNamespaceLoaded("stats4")) stats4::nobs   else nobs
  if(!missing(...)) {# several objects: produce data.frame
    lls <- lapply(list(object, ...), ll)
    vals <- sapply(lls, function(el) {
      no <- attr(el, "nobs")
      c(as.numeric(el), attr(el, "df"),
        if(is.null(no)) NA_integer_ else no)
    })
    val <- data.frame(df = vals[2L,], ll = vals[1L,], nobs = vals[3L,])
    nos <- na.omit(val$nobs)
    if (length(nos) && any(nos != nos[1L]))
      warning("models are not all fitted to the same number of observations")
    ## if any val$nobs = NA, try to get value via nobs().
    unknown <- is.na(val$nobs)
    if(any(unknown))
      val$nobs[unknown] <-
      sapply(list(object, ...)[unknown],
             function(x) tryCatch(Nobs(x), error = function(e) NA_real_))
    val <- data.frame(df = val$df, AICc = -2*val$ll + 2*val$df + 2*val$df*(val$df+1)/(val$nobs-val$df-1))
    row.names(val) <- as.character(match.call()[-1L])
    val
  } else {
    lls <- ll(object)
    nos <- attr(lls, "nobs")
    if (is.null(nos)) ## helps if has nobs() method, but logLik() gives no "nobs":
      nos <- tryCatch(Nobs(object), error = function(e) NA_real_)
    -2 * as.numeric(lls) + 2 * attr(lls, "df") + 2*attr(lls, "df")*(attr(lls, "df")+1)/(nos - attr(lls, "df") - 1)
  }
}

#' Plot method for tscopulafit class
#'
#' @param x an object of class \linkS4class{tscopulafit}.
#' @param plottype type of plot required.
#' @param bw logical variable specifying whether black-white options should be chosen.
#' @param lagmax maximum lag value for Kendall plots
#'
#' @return No return value, generates plot.
#' @export
#'
#' @examples
#' data <- sim(armacopula(list(ar = 0.5, ma = 0.4)), n = 1000)
#' fit <- fit(armacopula(list(ar = 0.5, ma = 0.4)), data)
#' plot(fit)
setMethod("plot", c(x = "tscopulafit", y = "missing"),
          function(x, plottype = "residual", bw = FALSE, lagmax = 30) {
            colchoice <- ifelse(bw, "gray50", "red")
            if ((plottype == "vtransform") & (is(x@tscopula, "vtscopula"))){
              plot(x@tscopula@Vtransform)
              U <- as.numeric(x@data)
              V <- vtrans(x@tscopula@Vtransform, U)
              points(strank(U), strank(V), col = colchoice)
              lines(sort(U), V[order(U)])
            }
            else if (plottype == "residual"){
              res <- resid(x)
              qqnorm(res)
              abline(0,1, col = colchoice)
            }
            else if (plottype == "kendall"){
              tauE <- glag(x, lagmax)
              k <- length(tauE)
              tauT <- kendall(x@tscopula, k)
              plot(1:k, tauE, type = "h", xlim = c(1, k),
                   ylim = range(tauE,tauT,0), xlab = "lag", ylab = "tau")
              abline(h=0)
              if (k >1)
                lines(1:k, tauT, col = colchoice)
              else
                points(1, tauT, col = colchoice)
            }
            else if (plottype == "glag"){
              ldata <- glag(x, lagmax, glagplot = TRUE)
              nplots <- length(ldata)
              lc <- ifelse(nplots > 4, 3, 2)
              lr <- ceiling(nplots / lc)
              oldpar <- par(no.readonly = TRUE)
              on.exit(par(oldpar))
              par(mfrow = c(lr, lc),
                  mar = c(2.1, 2.1, 1.5, 0.5),
                  oma = rep(2, 4),
                  pty = "s", cex = 0.5)
              for (i in 1:nplots)
                plot(ldata[[i]], main = paste("Lag ", i, sep = ""), asp = 1,
                     xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "")
            }
            else
              stop("Not a valid plot type")
})

#' @describeIn tscopulafit Residual method for tscopulafit class
#'
#' @param object an object of class \linkS4class{tscopulafit}.
#' @param trace extract trace instead of residuals.
#'
#' @export
#'
setMethod("resid", "tscopulafit",
          function(object, trace = FALSE) {
            copula <- object@tscopula
            data <- object@data
            if (is(copula, "vtscopula")){
              data <- vtrans(copula@Vtransform, data)
              copula <- copula@Vcopula
            }
            residfunc <- eval(parse(text=paste("resid_",is(copula)[1],sep="")))
            residfunc(copula, data, trace)
          })

#' Generalized lagging function
#'
#' @param x an object of class \linkS4class{tscopulafit}.
#' @param lagmax maximum value for lag.
#' @param glagplot logical value indicating generalized lag plot.
#'
#' @return If \code{glagplot} is \code{TRUE} a list of generalized lagged datasets
#' of maximum length 9 is returned to facilitate a generalized lagplot.
#' If \code{glagplot} is \code{FALSE} a vector of length \code{lagmax} containing
#' the Kendall rank correlations for the generalized lagged datasets is returned.
#' @export
#'
glag <- function(x, lagmax = 20, glagplot = FALSE) {
  data <- x@data
  copula <- x@tscopula
  if (is(copula, "vtscopula")){
    data <- vtrans(copula@Vtransform, data)
    copula <- copula@Vcopula
  }
  coptype <- is(copula)[1]
  lagfunc <- eval(parse(text=paste("glag_for_", coptype, sep="")))
  lagfunc(copula, data, lagmax, glagplot)
}

#' @describeIn tscopulafit Prediction method for tscopulafit class
#'
#' @param object an object of class \linkS4class{tscopulafit}.
#' @param x vector of arguments of prediction function.
#' @param type type of prediction function ("df" for density, "qf" for quantile function
#' or "dens" for density).
#'
#' @export
#'
setMethod("predict", c(object = "tscopulafit"), function(object, x, type = "df") {
  predict(object@tscopula, object@data, x, type = type)
})

#' Transform a fitted armacopula into a fitted dvinecopula or dvinecopula2 object
#'
#' @param object an object of class \linkS4class{tscopulafit} in which the copula is
#' of class \linkS4class{armacopula}.
#'
#' @return An object of class \linkS4class{tscopulafit} in which the copula is
#' a \linkS4class{dvinecopula} (for fitted AR copulas)
#' or class \linkS4class{dvinecopula2} (for fitted MA or ARMA copulas).
#' @export
armafit2dvine <- function(object){
  if (!(is(object, "tscopulafit")))
    stop("Not fitted armacopula")
  tscop <- arma2dvine(object@tscopula)
  new("tscopulafit", tscopula = tscop, data = object@data, fit = object@fit)
}
