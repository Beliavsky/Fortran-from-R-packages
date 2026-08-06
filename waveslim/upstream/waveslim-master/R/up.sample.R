#' Upsampling of a vector
#' 
#' Upsamples a given vector.
#' 
#' 
#' @param x vector of observations
#' @param f frequency of upsampling; e.g, 2, 4, etc.
#' @param y value to upsample with; e.g., NA, 0, etc.
#' @return A vector twice its length.
#' @author B. Whitcher
#' @references Any basic signal processing text.
#' @keywords ts
#' @export up.sample
up.sample <- function(x, f, y=NA) {
  n <- length(x)
  as.vector(rbind(x, matrix(rep(y, (f-1)*n), nrow=f-1)))
}
