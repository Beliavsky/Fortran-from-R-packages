#' Specify a Two-Dimensional Monte Carlo Model for Looped Evaluation
#'
#' Builds a \samp{mcmodelcut} object that can be sent to \code{\link{evalmccut}}.
#'
#' @rdname mccut
#' @param x for \code{mcmodelcut}: a call or an expression (if \code{is.expr=TRUE}) including
#'   three blocks (see Details and Examples for the required structure). For \code{print.mccut}
#'   and \code{plot.mccut}: an \samp{mccut} object.
#' @param is.expr \code{FALSE} to send a call, \code{TRUE} to send an expression
#'   (see \code{\link{mcmodel}} examples).
#' @export
mcmodelcut <- function(x, is.expr=FALSE)
{
  if(!is.expr) x <- as.expression(substitute(x))
  if(!is.expression(x)) stop("x can not be evaluate as an expression")

  nbexpr <- length(x[[1]])
  if(nbexpr != 4) stop("The expression should include three blocks.")

  cat("The following expression will be evaluated only once :\n")
  cat(deparse(x[[1]][[2]]),sep="\n")

  last <- x[[1]][[3]][[length(x[[1]][[3]])]][[3]]
  lastcall1 <- substr(deparse(last,nlines=1), 1, 3)
  if (lastcall1 != "mc(") warning("The last call should be 'mymc <- mc(...)'")

  nom <- as.character(x[[1]][[3]][[length(x[[1]][[3]])]][[2]])                  # name of the mc object in the third block
  cat("The mc object is named: ",nom,"\n")

  class(x) <- "mcmodelcut"
  return(invisible(x))
}
