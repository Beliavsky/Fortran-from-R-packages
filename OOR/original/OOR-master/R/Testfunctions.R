##' @name Test functions
##' @aliases sin1
##' @aliases difficult
##' @description Several test functions of varying complexity are available. They are defined on [0,1].
##' @title Test functions of \code{x}
##' @param x vector specifying the location where the function is to be evaluated.
##' @param rho1,rho2,tmax additional parameters for double_sine.
##' @details These test functions are translated from the Matlab and Python codes in the references.
##' @references
##' M. Valko, A. Carpentier and R. Munos (2013), Stochastic Simultaneous Optimistic Optimization,
##' \emph{ICML}, 19-27 \url{https://inria.hal.science/hal-00789606}. Matlab code: \url{https://team.inria.fr/sequel/software/StoSOO/}. \cr \cr
##' J.-B. Grill, M. Valko and R. Munos (2015), Black-box optimization of noisy functions with unknown smoothness,
##'   \emph{NIPS}, 667-675 \url{https://inria.hal.science/hal-01222915}. Python code: \url{https://team.inria.fr/sequel/software/POO/}. \cr \cr
##' @examples
##' par(mfrow = c(2,3))
##'
##' curve(guirland, n = 501)
##' curve(sin1)
##' curve(difficult, xlim = c(1e-8, 1), n = 1001)
##' xgrid <- seq(0, 1, length.out = 500)
##' plot(xgrid, sapply(xgrid, difficult2), type = 'l', ylab = "difficult2(x)")
##' plot(xgrid, sapply(xgrid, double_sine), type = 'l', ylab = "double_sine(x) (default)")
##' double_sine2 <- function(x) double_sine(x, rho1 = 0.8, rho2 = 0.3)
##' plot(xgrid, sapply(xgrid, double_sine2), type = 'l', ylab = "double_sine(x) (modified)")
##'
##' par(mfrow = c(1,1))




##' @rdname  Testfunctions
##' @export
guirland <- function(x){
  return(4*x*(1-x)*(0.75+0.25*(1-sqrt(abs(sin(60*x))))))
}


##' @rdname  Testfunctions
##' @export
sin1 <- function(x){
  return((sin(13 * x) * sin(27 * x) / 2.0 + 0.5))
}


##' @rdname  Testfunctions
##' @export
difficult <- function(x){
  return(1-sqrt(x) + (-x*x +sqrt(x) )*(sin(1/(x*x*x))+1)/2)
}


##' @rdname Testfunctions
##' @export
difficult2 <- function(x){
  tmp <- log2(abs(x - 0.5))
  if(tmp %% 1 <= 0.5){
    tmp <- 1
  }else{
    tmp <- 0
  }
  return(tmp * (sqrt(abs(x - 0.5)) - (x - 0.5)^2) - sqrt(abs(x - 0.5)))
}


##' @rdname Testfunctions
##' @export
double_sine <- function(x, rho1 = 0.3, rho2 = 0.8, tmax = 0.5){
  u <- 2 * abs(x - tmax)
  if(u == 0) return(0)
  enveloppe_width <- u^(-log2(rho2)) - u^(-log2(rho1))
  return(mysin2(log2(u)/2) * enveloppe_width - u^(-log2(rho2)))
}

mysin2 <- function(x) return((sin(2*x*pi) + 1)/2)

