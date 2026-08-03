# File: GenSA.R
# Author: Sylvain Gubian
# Aim: Function for General Simulated Annealing

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

###############################################################################
#' Generalized Simulated Annealing Function
#' @description
#' This function searches for global minimum of a very complex non-linear
#' objective function with a very large number of optima.
#' @details
#' The default values of the control components are set for a complex
#' optimization problem.
#' For usual optimization problem with medium complexity, GenSA can find a
#' reasonable solution quickly sot he user is recommended to let GenSA stop
#' earlier by setting \code{threshold.stop}. If \code{threshold.stop} is the
#' expected function value, or by setting \code{max.time}. If the user just
#' want to run GenSA for \code{max.time} seconds, or by setting \code{max.call}.
#' If the user just want to run GenSA within \code{max.call} function calls.
#' Please refer to the examples below. For very complex optimization problems,
#' the user is recommended to increase \code{maxit} and \code{temp}.
#' @param par Vector. Initial values for the components to be optimized.
#' Default is \code{NULL}, in which case, default values will be generated
#' automatically.
#' @param fn A function to be minimized, with first argument the vector of
#' parameters over which minimization is to take place. It should return
#' a scalar result.
#' @param lower Vector with length of \code{par}. Lower bounds for components.
#' @param upper Vector with length of \code{par}. Upper bounds for components.
#' @param ... allows the user to pass additional arguments to the function
#' \code{fn}.
#' @param control The argument is a list that can be used to control the
#' behavior of the algorithm
#'    \describe{
#'        \item{\code{maxit}}{
#'            Integer. Maximum number of iterations of the algorithm.
#'        }
#'        \item{\code{threshold.stop}}{
#'            Numeric. The program will stop when the expected objective
#'            unction value \code{threshold.stop} is reached. Default value
#'            is \code{NULL}
#'        }
#'        \item{\code{nb.stop.improvement}}{
#'            Integer. The program will stop when there is no any improvement
#'            in \code{nb.stop.improvement} steps.
#'        }
#'        \item{\code{smooth}}{
#'              Logical.\code{TRUE} when the objective function is smooth, or
#'              differentiable almost everywhere in the region of \code{par},
#'              \code{FALSE} otherwise. Default value is \code{TRUE}.
#'          }
#'          \item{\code{max.call}}{
#'              Integer. Maximum number of call of the objective function.
#'              Default is set to 1e7.
#'          }
#'          \item{\code{max.time}}{
#'              Numeric. Maximum running time in seconds.
#'          }
#'          \item{\code{temperature}}{
#'              Numeric. Initial value for temperature.
#'          }
#'          \item{\code{visiting.param}}{
#'              Numeric. Parameter for visiting distribution.
#'          }
#'          \item{\code{acceptance.param}}{
#'              Numeric. Parameter for acceptance distribution.
#'          }
#'          \item{\code{verbose}}{
#'              Logical. \code{TRUE} means that messages from the algorithm are
#'              shown. Default is \code{FALSE}.
#'          }
#'          \item{\code{simple.function}}{
#'              Logical. \code{FALSE} means that the objective function has only
#'              a few local minima. Default is \code{FALSE} which means that the
#'              objective function is complicated with many local minima.
#'          }
#'          \item{\code{trace.mat}}{
#'              Logical. Default is \code{TRUE} which means that the trace
#'              matrix will be available in the returned value of \code{GenSA}
#'              call.
#'          }
#'          \item{\code{seed}}{
#'              Integer. Negative integer value that can be set to
#'              initialize the internal random generator.
#'          }
#'          }
#'@returns The returned value is a list with the following fields:
#'          \describe{
#'              \item{value}{
#'                  Numeric. The value of \code{fn} corresponding to \code{par}.
#'              }
#'              \item{par}{
#'                  Vector. The best set of parameters found.
#'              }
#'              \item{trace.mat}{
#'                  A matrix which contains the history of the algorithm.
#'                  (By columns: Step number, temperature,
#'                  current objective function value, current minimal objective
#'                  function value).
#'              }
#'              \item{counts}{
#'                  Integer. Total number of calls of the objective function.
#'              }
#'          }
#' @references{
#'  Xiang Y, Gubian S, Martin F (2017). "Generalized Simulated Annealing."
#'  IntechOpen, Computational Optimization in Engineering, Chapter 2.
#'
#'  Xiang Y, Gubian S, Suomela B, Hoeng (2013). "Generalized Simulated Annealing
#'  for Efficient Global Optimization: the GenSA Package for R". The R Journal
#'  Volume 5/1, June 2013.
#'
#'  Xiang Y, Sun DY, Gong XG (2000). "Generalized Simulated Annealing Studies on
#'  Structures and Properties of Nin (n=2-55) Clusters." Journal of Physical
#'  Chemistry A, 104, 2746--2751.
#'
#'  Xiang Y, Gong XG (2000a). "Efficiency of Generalized Simulated Annealing."
#'  PHYSICAL REVIEW E, 62, 4473.
#'
#'  Xiang Y, Sun DY, Fan W, Gong XG (1997). "Generalized Simulated Annealing
#'  Algorithm and Its Application to the Thomson Model." Physics Letters A, 233,
#'  216--220.
#'
#'  Tsallis C, Stariolo DA (1996). "Generalized Simulated Annealing."
#'  Physica A, 233, 395--406.
#'
#'  Tsallis C (1988). "Possible generalization of Boltzmann-Gibbs statistics."
#'  Journal of Statistical Physics, 52, 479--487.
#'}
#'@examples
#'  library(GenSA)
#'  # Try Rastrgin function (The objective function value for global minimum
#'  # is 0 with all components of par are 0.)
#'  Rastrigin <- function(x) {
#'      sum(x^2 - 10 * cos(2 * pi  * x)) + 10 * length(x)
#'  }
#'  # Perform the search on a 30 dimensions rastrigin function. Rastrigin
#'  # function with dimension 30 is known as the most
#'  # difficult optimization problem according to "Yao X, Liu Y, Lin G (1999).
#'  # \Evolutionary Programming Made Faster."
#'  # IEEE Transactions on Evolutionary Computation, 3(2), 82-102.

#'  # GenSA will stop after finding the targeted function value 0 with
#'  # absolute tolerance 1e-13
#'  set.seed(1234) # The user can use any seed.
#'  dimension <- 30
#'  global.min <- 0
#'  tol <- 1e-13
#'  lower <- rep(-5.12, dimension)
#'  upper <- rep(5.12, dimension)
#'  out <- GenSA(lower = lower, upper = upper, fn = Rastrigin,
#'               control=list(threshold.stop=global.min+tol,verbose=TRUE))
#'  out[c("value","par","counts")]
#'
#'  # GenSA will stop after running for about 2 seconds
#'  # Note: The time for solving this problem by GenSA may vary
#'  # depending on the computer used.
#'  set.seed(1234) # The user can use any seed.
#'  dimension <- 30
#'  global.min <- 0
#'  tol <- 1e-13
#'  lower <- rep(-5.12, dimension)
#'  upper <- rep(5.12, dimension)
#'  out <- GenSA(lower = lower, upper = upper, fn = Rastrigin,
#'               control=list(max.time=2))
#'  out[c("value","par","counts")]
#'

#'@author {
#'    Yang Xiang, Sylvain Gubian, Brian Suomela, Julia Hoeng, PMP SA.
#' .  (Y.Xiang and S.Gubian are equal contributors)
#'}

#'@export
#'@importFrom stats constrOptim runif
#'@importFrom utils write.table
#' @useDynLib GenSA, createInstance, releaseInstance, execute, getREnergy,
#' getRXMiniVector, getRTraceMat, getRNbFuncCall,getRTraceMatSize

GenSA <- function(par = NULL, fn, lower, upper, control = list(), ...) {
    # Do some checks
    jc <- NULL
    if (!is.function(fn) || is.null(fn)) {
        stop("'fn' has to be a R function.")
    }

    # Create an environment for sharing data between R and C
    gensa_env <- new.env(hash = TRUE, parent = emptyenv())

    fn1 <- function(par) {
        ret <- fn(par, ...)
        return(ret)
    }


    if (!is.null(jc)) {
        if (!is.function(jc)) {
            stop("'jc' has to be a R function.")
        }
        jc1 <-  function(par, ...) {
            return(jc(par, ...))
        }
    } else {
        jc1 <- NULL
    }

    LSE <- function(theta, ui, ci, mu, xlow, xhigh, count) {
        assign("xlow", xlow, envir = gensa_env)
        assign("xhigh", xhigh, envir = gensa_env)
        assign("count", count, envir = gensa_env)

        res <- constrOptim(theta = theta, f = fn2, ui = ui, ci = ci, mu = mu,
            grad = NULL, outer.eps = 1e-06)
        counts <- get("count", gensa_env)
        ret <- list(value = res$value, convergence = res$convergence,
            par = res$par, counts = as.integer(counts))
        return(ret)
    }

    fn2 <- function(par, ...) {
        if (!is.null(jc1)) {
            in_constraint <- jc1(par, ...)
            if (!in_constraint) {
                return(1e10)
            }
        } else {
            penalty <- 0
            xlow <- get("xlow", gensa_env)
            xhigh <- get("xhigh", gensa_env)
            count <- get("count", gensa_env)
            counts <- count
            for (i in (length(lower))) {
                if (par[i] >= xlow[i] && par[i] <= xhigh[i]) {
                    delta_energy <- 0
                } else {
                    if (par[i] < xlow[i]) {
                        delta_energy <-  abs(par[i] - xlow[i]) * 1e11
                    }
                    if (par[i] > xhigh[i]) {
                        delta_energy <- abs(par[i] - xhigh[i]) * 1e11
                    }
                }
                penalty <- penalty + delta_energy
            }
            if (penalty > 1.e-10) {
                to_return <- penalty+1.e10
                return(to_return)
            } else {
                to_return <- fn1(par, ...)
                counts <- counts + 1
                to_return <- to_return + penalty
                if (is.nan(to_return)) {
                    to_return <- count * 1e5 + 1e10
                    assign("count", counts, envir = gensa_env)
                    return(to_return)
                }  else {
                    assign("count", counts, envir = gensa_env)
                    return(to_return)
                }
            }
        }
    }

    assign("LSE", LSE, envir = gensa_env)

    con <- list(
                maxit = 5000,
                threshold.stop = NULL,
                temperature = 5230,
                visiting.param = 2.62,
                acceptance.param = -5.0,
                max.time = NULL,
                nb.stop.improvement = -1,
                smooth = TRUE,
                max.call = 10000000,
                simple.function = FALSE,
                trace.fn = NULL,
                verbose = FALSE,
                trace.mat = TRUE,
                seed = -100377
                )
    con$high.dim <- TRUE
    con$markov.length <- 2 * length(lower)
    con$tem.restart <- .1

    # Perform some checks before callinc C code
    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC]))
        warning("unknown names in control: ", paste(noNms, collapse = ", "))


    if (!exists("par") && (length(lower) == 0 || length(upper) == 0)) {
        stop("There is no par or no lower/upper bounds defined")
    }

    if (length(lower) != length(upper)) {
        stop("Lower and upper bounds vector do not have the same length")
    }

    cmp <- unique(lower < upper)
    if (length(cmp) != 1 || !cmp) {
        stop("Lower and upper bounds are not consistent (lower >= upper)")
    }

    if (!is.null(par)) {

        if (length(lower) == 0 || length(lower) != length(par)) {
            stop("Lower bounds vector size does not match with par size")
        }
        if (length(upper) == 0 || length(upper) != length(par)) {
            stop("Upper bounds vector size does not match with par size")
        }
        if (any(is.na(par)) || any(is.nan(par)) || any(is.infinite(par))) {
            stop("par contains NA, NAN or Inf")
        }
    } else {
        if (con$verbose) {
            cat("Initializing par with random data inside bounds\n")
        }
        par <- vector()
        #initialize par with random values in the bounds
        par <- lower + runif(length(lower)) * (upper - lower)
    }

    ret <- list()
    # Create instance of the GenSACaller
    instance <- .Call("createInstance", PACKAGE = "GenSA")
    if (is.null(instance)) {
        stop("Can not create GenSACaller instance!")
    }

    # Call execute on the instance
    res <- .Call("execute", par, lower, upper, fn1, jc1, con, gensa_env,
        instance, PACKAGE = "GenSA")
    if (is.null(res)) {
        stop("Can not call execute function on instance")
    }

    # Get the results in a list
    res <- .Call("getREnergy", instance, PACKAGE = "GenSA")
    if (is.null(res)) {
        message("Can not get minimum function value")
    } else {
        ret$value <- res
    }

    res <- .Call("getRXMiniVector", instance, PACKAGE="GenSA")
    if (is.null(res)) {
        message("Can not get calculated par values")
    } else {
        ret$par <- res
    }

    if (con$trace.mat) {
        nr <- .Call("getRTraceMatSize", instance, PACKAGE = "GenSA")
        if (nr > 0) {
            ret$trace.mat <- matrix(NA, nr, 4)
            ret$trace.mat[, 1] <- as.integer(.Call("getRTraceMat", instance,
                "nSteps", PACKAGE = "GenSA"))
            ret$trace.mat[, 2] <- as.numeric(.Call("getRTraceMat", instance,
                "temperature", PACKAGE = "GenSA"))
            ret$trace.mat[,3] <- as.numeric(.Call("getRTraceMat", instance,
                "currentEnergy", PACKAGE = "GenSA"))
            ret$trace.mat[,4] <- as.numeric(.Call("getRTraceMat", instance,
                "minEnergy", PACKAGE = "GenSA"))
            colnames(ret$trace.mat) <- c("nb.steps", "temperature",
                "function.value", "current.minimum")
        }
        if (!is.null(con$trace.fn)) {
            if (con$verbose) {
                cat(paste("Writing trace.mat data into file:",
                    con$trace.fn, "\n"))
            }
            write.table(ret$trace.mat, file = con$trace.fn, col.names = TRUE,
                row.names = FALSE)
            ret$trace.mat <- paste("trace.mat is written in file:",
                con$trace.fn)
        }
    }


    res <- .Call("getRNbFuncCall", instance, PACKAGE = "GenSA")
    if (is.null(res)) {
        message("Can not get number of function calls")
    } else {
        ret$counts <- res
    }
    .Call("releaseInstance", instance, PACKAGE = "GenSA")
    ret
}
