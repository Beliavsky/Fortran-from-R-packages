#' Build mcnode Objects Without Dimension Control
#'
#' \samp{mcdatanocontrol} is a fast but unchecked version of \code{\link{mcdata}} which
#' forces the dimension of data to be \samp{(nsv x nsu x nvariates)} and sets the attributes
#' and class without any validation. Useful when your model is tested and speed is important.
#'
#' @rdname mcnode
#' @export
mcdatanocontrol <- function(data, type=c("V","U","VU","0"), nsv=ndvar(), nsu=ndunc(), nvariates=1, outm="each")
{
  dim(data) <- NULL
  data[1:(nsv*nsu*nvariates)] <- data
  dim(data) <- c(nsv,nsu,nvariates)
  class(data) <- "mcnode"
  attr(data,which="type") <- type
  attr(data,which="outm") <- outm
  }
