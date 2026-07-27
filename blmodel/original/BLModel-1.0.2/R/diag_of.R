#' Extracts the diagonal from  matrix x.
#'
#'@param x g
#'
#'

.diag_of <- function (x)
{
  if ((m <- min(dim(x))) == 0)
    return(numeric(0))
  y <- c(x)[1 + 0:(m - 1) * (dim(x)[1] + 1)]
  nms <- dimnames(x)
  if (is.list(nms) && !any(sapply(nms, is.null)) && identical((nm <- nms[[1]][1:m]),
                                                              nms[[2]][1:m]))
    names(y) <- nm
  return(y)
}
