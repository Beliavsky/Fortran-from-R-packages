#' Evaluates a Monte Carlo Model
#'
#' Evaluates a \code{\link{mcmodel}} object (or a valid expression) using a specified number
#' of simulations and with (or without) a specified seed.
#'
#' The model is evaluated. Intermediate variables used to build the \samp{mc} object are not
#' stored.
#'
#' @note The seed is set at the beginning of the evaluation. Thus, the complete similarity
#' of two evaluations with the same seed is not certain, depending on the structure of the model.
#'
#' @param expr a model of class \code{\link{mcmodel}} or a valid expression.
#' @param nsv the number of simulations in the variability dimension.
#' @param nsu the number of simulations in the uncertainty dimension.
#' @param seed the random seed. If \samp{NULL} the seed is unchanged.
#' @return The result of the evaluation — should be a \samp{mc} object.
#' @seealso \code{\link{mcmodel}}, \code{\link{evalmccut}} to evaluate high-dimension models.
#' @keywords methods
#' @examples
#' data(ec)
#' ec$modEC1
#' evalmcmod(ec$modEC1, nsv=100, nsu=100, seed=666)
#' @export
evalmcmod <- function(expr, nsv = ndvar(), nsu = ndunc(), seed = NULL)
{
  if(!is.null(seed)) set.seed(seed)
  Oldv <- ndvar()
  Oldu <- ndunc()
  ndvar(nsv)
  ndunc(nsu)
  x <- try({
	  x <- eval(expr)
	  if(!is.mc(x)) stop("expr does not lead to a mc")
    x}, silent=TRUE)

  ndvar(Oldv)
  ndunc(Oldu)
  if(inherits(x,"try-error")) stop(x,call. = FALSE)
  return(x)
}
