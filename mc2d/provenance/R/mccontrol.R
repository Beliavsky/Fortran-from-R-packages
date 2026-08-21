.pkgglobalenv <- new.env(parent=emptyenv())

#' Sets or Gets the Default Number of Simulations
#'
#' Sets or retrieves the default number of simulations in the variability and
#' uncertainty dimensions.
#'
#' \samp{ndvar()} gets and \samp{ndvar(n)} sets the default number of simulations in
#' the 1D simulations or the number of simulations in the variability dimension in
#' 2D simulations.
#'
#' \samp{ndunc()} gets and \samp{ndunc(n)} sets the number of simulations in the
#' uncertainty dimension in 2D simulations.
#'
#' \samp{n} is rounded to its ceiling value.
#'
#' The default values when loaded are 1001 for \samp{ndvar} and 101 for \samp{ndunc}.
#'
#' @param n number of simulations.
#' @return The current value, AFTER modification if \samp{n} is present (unlike \code{options}).
#' @seealso \code{\link{mcstoc}}, \code{\link{mcdata}}
#' @keywords misc
#' @name mc.control
#' @examples
#' (oldvar <- ndvar())
#' (oldunc <- ndunc())
#' mcstoc(runif, type = "VU")
#' ndvar(12)
#' ndunc(21)
#' mcstoc(runif, type = "VU")
#' ndvar(oldvar)
#' ndunc(oldunc)
#' @export
ndvar <- function(n)
{
  if(!exists("mc.control", envir=.pkgglobalenv))
    assign("mc.control",list(nsv=1001,nsu=101),envir=.pkgglobalenv )
  x <- get("mc.control", envir=.pkgglobalenv)
   if(!is.list(x) || is.null(x$nsv) || is.null(x$nsu))
    assign("mc.control",list(nsv=1001,nsu=101),envir=.pkgglobalenv )
  if(!missing(n)){
    if (n > 0) x$nsv <- ceiling(n)
        else stop("Invalid n")
    assign("mc.control",x, envir=.pkgglobalenv)}
  return(x$nsv)}

#' @rdname mc.control
#' @export
ndunc <- function(n)
{
  if(!exists("mc.control", envir=.pkgglobalenv))
    assign("mc.control",list(nsv=1001,nsu=101),envir=.pkgglobalenv )
  x <- get("mc.control", envir=.pkgglobalenv)
   if(!is.list(x) || is.null(x$nsv) || is.null(x$nsu))
    assign("mc.control",list(nsv=1001,nsu=101),envir=.pkgglobalenv )
  if(!missing(n)){
    if (n > 0) x$nsu <- ceiling(n)
        else stop("Invalid n")
    assign("mc.control",x, envir=.pkgglobalenv)}
  return(x$nsu)}
