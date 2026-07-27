#'Makes a diagonal matrix with values from x
#'
#'
#'Function make_diag() takes a vector x and returns a matrix with diagonal consists the values of x.
#'
#'@param x vector
#'@return diagonal matrix y
#'@keywords internal 


.make_diag <- function (x){
  if(class(x)!= "numeric"){
    stop(paste("Wrong input"))
  }
  n <- length(x)
  y <- array(0, c(n, n))
  if (n > 0)
    y[1 + 0:(n - 1) * (n + 1)] <- x
  return(y)
}
