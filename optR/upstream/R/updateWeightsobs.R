# Function to transform the x and y for weighted regression
updateWeightsobs<-function(x, y, w){
  if (is.null(n <- nrow(x))) 
    stop("'x' must be a matrix")
  if (n == 0) 
    stop("0 (non-NA) cases")
  if (is.matrix(y) && ny == 1L) 
    y <- drop(y)
  
  ny <- NCOL(y)
  zero.weights <- any(w == 0)
  if (zero.weights) {
    ok <- w != 0
    nok <- !ok
    w <- w[ok]
    x0 <- x[!ok, , drop = FALSE]
    x <- x[ok, , drop = FALSE]
    n <- nrow(x)
    
    y0 <- if (ny > 1L) 
      y[!ok, , drop = FALSE]
    else y[!ok]
    
    y <- if (ny > 1L)
      y[ok, , drop = FALSE]
    else y[ok]
  }
  
  # Update the weights
  wts<-sqrt(w)
  x<-x*wts
  y<-y*wts
  
  return(list("x"=x, "y"=y, "weights"=w))
}