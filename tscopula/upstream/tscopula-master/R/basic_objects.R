#' Time series copula processes
#'
#' Class of objects for time series copula processes.
#'
#'
#' @export
#'
#' @import stats
#' @import stats4
#' @import methods
#' @import graphics
#' @import utils
#' @import zoo
#' @import xts
#'
setClass("tscopula", contains = c("VIRTUAL"))

#' Generic for simulating time series copula models
#'
#' Methods are available for objects of class \linkS4class{swncopula},
#' \linkS4class{armacopula},
#' \linkS4class{dvinecopula}, \linkS4class{dvinecopula2},
#' \linkS4class{margin} and \linkS4class{tscm}.
#'
#' @param object an object of the model class.
#' @param ... further arguments to be passed to the simulation.
#'
#' @return A simulated realization from the time series model.
#' @export
#'
#'
setGeneric("sim", function(object, ...) {
  standardGeneric("sim")
})

#' Strict white noise copula process
#'
#' @export
#'
setClass("swncopula", contains = "tscopula")

#' Constructor function for strict white noise copula process
#'
#' @return Object of class \linkS4class{swncopula}.
#' @export
#'
#' @examples
#' swncopula()
swncopula <- function() {
  new("swncopula")
}

#' @describeIn swncopula Simulation method for strict white noise copula
#'
#' @param object an object of the class.
#' @param n numeric value for length of simulated realisation.
#'
#' @export
#'
#' @examples
#' sim(swncopula())
setMethod("sim", "swncopula", function(object, n = 1000) {
  runif(n)
})

#' @describeIn swncopula Coef method for strict white noise copula
#'
#' @param object an object of the class.
#'
#' @export
#'
setMethod("coef", "swncopula", function(object) {
  NULL
})

#' @describeIn swncopula Show method for strict white noise copula
#'
#' @param object an object of class \linkS4class{swncopula}.
#'
#' @export
#'
setMethod("show", "swncopula", function(object) {
  cat("SWN \n")
})

#' Generic for Kendall correlations
#'
#' Methods are available for objects of class \linkS4class{armacopula},
#' \linkS4class{dvinecopula}, \linkS4class{dvinecopula2} and \linkS4class{vtscopula}.
#'
#' @param object an object of the model class.
#' @param ... further arguments to be passed to Kendall calculation.
#'
#' @return A vector of Kendall correlations.
#' @export
#'
#'
setGeneric("kendall", function(object, ...) {
  standardGeneric("kendall")
})



