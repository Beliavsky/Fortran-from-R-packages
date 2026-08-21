#' Specify a Monte Carlo Model
#'
#' Specifies a \samp{mcmodel} without evaluating it, for subsequent evaluation using
#' \code{\link{evalmcmod}}.
#'
#' The model should be put between \samp{\{} and \samp{\}}. The last line should be of the
#' form \samp{mc(...)}. Any reference to the number of simulations in the variability
#' dimension should use \samp{ndvar()} or (preferred) \samp{nsv}. Any reference to the
#' number of simulations in the uncertainty dimension should use \samp{ndunc()} or (preferred)
#' \samp{nsu}.
#'
#' @param x an \R call or an expression.
#' @param is.expr \samp{FALSE} to send a call, \samp{TRUE} to send an expression.
#' @return an \R expression of class \samp{mcmodel}.
#' @seealso \code{\link{expression}}, \code{\link{evalmcmod}} to evaluate the model,
#'   \code{\link{mcmodelcut}} to evaluate in a loop.
#' @keywords methods
#' @examples
#' modEC1 <- mcmodel({
#'   conc <- mcdata(10, "0")
#'   cook <- mcstoc(rempiricalD, values=c(0, 1/5, 1/50), prob=c(0.027, 0.373, 0.600))
#'   serving <- mcstoc(rgamma, shape=3.93, rate=0.0806)
#'   expo <- conc * cook * serving
#'   dose <- mcstoc(rpois, lambda=expo)
#'   risk <- 1 - (1 - 0.001)^dose
#'   mc(conc, cook, serving, expo, dose, risk)
#' })
#' evalmcmod(modEC1, nsv=100, nsu=100)
#' @export
mcmodel <- function(x, is.expr=FALSE)
{
  if(!is.expr) x <- as.expression(substitute(x))
  if(!is.expression(x)) stop("x can not be evaluate as an expression")

    last <- x[[1]][length(x[[1]])]
    lastcall1 <- substr(deparse(last, nlines=1), 1, 3)
    if (lastcall1 != "mc(") warning("The last call should be 'mc(...)'")
    class(x) <- "mcmodel"
    return(x)
}
