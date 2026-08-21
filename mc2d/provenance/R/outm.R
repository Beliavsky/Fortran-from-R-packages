#' Changes the Output of Nodes
#'
#' Changes the \samp{outm} attribute of an \samp{mcnode} or a node of an \samp{mc} object.
#'
#' @param x a \samp{mcnode} or a \samp{mc} object.
#' @param value the output of the \samp{mcnode} for multivariate nodes. May be \samp{"each"}
#'   (default) if output should be provided for each variate independently, \samp{"none"} for
#'   no output, or a vector of function names (as character strings) applied on the variates
#'   dimension before any output (e.g. \samp{"mean"}, \samp{"median"}, \samp{c("min","max")}).
#'   The function should return one value per vector.
#' @param which.node which node should be changed in a \samp{mc} object.
#' @return \samp{x} with a modified \samp{outm} attribute.
#' @keywords misc
#' @examples
#' data(total)
#' total$xVUM2
#' ## since outm = NULL
#' summary(total$xVUM2)
#' x <- outm(total$xVUM2, c("min"))
#' summary(x)
#' @export
outm <- function(x,value="each",which.node=1)
{
  if(is.character(value)  && (value == "none" || value == "each" | all(sapply(value,exists,mode="function")))){

  if(inherits(x,"mcnode"))  attr(x,which="outm") <- value
      else if(inherits(x,"mc")) {
        if(is.character(which.node)) which.node <- which(names(x) %in% which.node)
        attr(x[which.node],which="outm") <- value}
    }
  else stop("value should be 'none','each' or the name a valid function")
  return(x)
}
