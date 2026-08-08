#' This package implements optimistic optimization methods [1,2,3] for global optimization of deterministic or stochastic functions.
#' The algorithms feature guarantees of the convergence to a global optimum.
#' They require minimal assumptions on the (only local) smoothness, where the smoothness parameter does not need to be known.
#' They are expected to be useful for the most difficult functions when we have no information on smoothness and the gradients are unknown or do not exist.
#' Due to the weak assumptions, however, they can be mostly effective only in small dimensions, for example, for hyperparameter tuning [4].
#' @title Package OOR
#' @docType package
#' @name OOR
#' @aliases OOR-package
#' @references
#' [1] R. Munos (2011), Optimistic optimization of deterministic functions without the knowledge of its smoothness,
#' \emph{NIPS}, 783-791. \cr \cr
#' [2] M. Valko, A. Carpentier and R. Munos (2013), Stochastic Simultaneous Optimistic Optimization,
#' \emph{ICML}, 19-27 \url{https://inria.hal.science/hal-00789606}. \cr \cr
#' [3] J.-B. Grill, M. Valko and R. Munos (2015), Black-box optimization of noisy functions with unknown smoothness,
#'   \emph{NIPS}, 667-675 \url{https://inria.hal.science/hal-01222915}. \cr \cr
#' [4] S. Samothrakis, D. Perz, S. Lucas (2013), Training gradient boosting machines using curve-fitting and information-theoretic features for causal direction detection,
#' \emph{NIPS Workshop on Causality}.
#' @details
#' Important functions: \cr
#' \code{\link[OOR]{StoSOO}} \cr
#' \code{\link[OOR]{POO}} \cr
#' @note
#' This package is based on the Matlab and Python implementations from the corresponding publications,
#' available from the following webpage: \url{https://team.inria.fr/sequel/software/}.
#'
#' @examples
#' #------------------------------------------------------------
#' # Example 1 : Deterministic optimization with SOO
#' #------------------------------------------------------------
#' ## Define objective
#' fun1 <- function(x) return(-guirland(x))
#'
#' ## Optimization
#' Sol1 <- StoSOO(par = NA, fn = fun1, nb_iter = 1000, control = list(type = "det", verbose = 1))
#'
#' ## Display objective function and solution fund
#' curve(fun1, n = 1001)
#' abline(v = Sol1$par, col = 'red')
#'
#' #------------------------------------------------------------
#' # Example 2 : Stochastic optimization with StoSOO
#' #------------------------------------------------------------
#' set.seed(42)
#'
#' ## 2-dimensional noisy objective function, defined on [0, pi/4]^2
#' fun2 <- function(x){return(-sin1(x[1]) * sin1(1 - x[2]) + runif(1, min = -0.05, max = 0.05))}
#'
#' ## Optimizing
#' Sol2 <- StoSOO(par = rep(NA, 2), fn = fun2, upper = rep(pi/4, 2), nb_iter = 1000)
#'
#' ## Display solution
#' xgrid <- seq(0, pi/4, length.out = 101)
#' Xgrid <- expand.grid(xgrid, xgrid)
#' ref <- apply(Xgrid, 1, function(x){(-sin1(x[1]) * sin1(1 - x[2]))})
#' filled.contour(xgrid, xgrid, matrix(ref, 101), color.palette  = terrain.colors,
#' plot.axes = {axis(1); axis(2); points(Xgrid[which.min(ref),, drop = FALSE], pch = 21);
#'              points(Sol2$par[1], Sol2$par[2], pch = 13)})
#'
#' \dontrun{
#' #------------------------------------------------------------
#' # Example 3 : Stochastic optimization with POO
#' #------------------------------------------------------------
#' set.seed(10)
#' noise.level <- 0.05
#'
#' ## Define and display objective
#' fun3 <- function(x){return(double_sine(x) + runif(1, min = -noise.level, max = noise.level))}
#' xgrid <- seq(0, 1, length.out = 1000)
#' plot(xgrid, sapply(xgrid, double_sine), type = 'l', ylab = "double_sine(x)", xlab = 'x')
#'
#' ## Maximization
#' Sol3 <- POO(fun3, horizon = 1000, noise.level = noise.level)
#'
#' ## Display result
#' abline(v = Sol3$par)
#' }
NULL
