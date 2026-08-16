#' Draw the CDF for a discrete random variable
#'
#' @author Freddy Hernandez, \email{fhernanb@unal.edu.co}
#'
#' @param x vector with the values of the random variable \eqn{X}.
#' @param fx vector with the probabilities of \eqn{X}.
#' @param col color for the line.
#' @param lwd line width.
#' @param ... further arguments and graphical parameters.
#'
#' @return
#' A plot with the cumulative distribution function.
#'
#' @example  examples/examples_plot_cdf.R
#' @importFrom stats stepfun
#' @importFrom graphics grid
#' @export
#'
plot_discrete_cdf <- function(x, fx,
                              col="blue", lwd=3, ...) {
  Fx <- cumsum(fx)
  F <- stepfun(x=x, y=c(0, Fx), right=TRUE)

  # Para dibujar la funcion F(x)
  plot(F, verticals=FALSE,
       lwd=lwd, col=col, las=1,
       xlab="X", ylab="F(X=x)", ...)

  grid()   # Para incluir una rejilla
}
