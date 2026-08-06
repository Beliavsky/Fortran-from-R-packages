tv <- function (speed, location, xtv = NULL, n = NULL, opt = 0, order.g = NULL, 
                as.zoo = TRUE, verbose = FALSE)
{
  if( is.null(order.g) ){ 
    order.g <- length(location) 
  }
  if (length(location) != order.g) {
    stop("Check the number of location parameters.")
  }
  if (is.null(xtv)) { 
    if (!is.null(n)) { 
      xtv <- seq(1:n)/n
    }
  }
  xtv <- as.numeric(xtv)
  if (is.null(n)) n <- length(xtv)
  if (opt == 0) speedadj <- speed
  if (opt == 1) speedadj <- speed/sd(xtv)
  if (opt == 2) speedadj <- exp(speed)
  mat <- matrix(rep(location, n), n, order.g, byrow = TRUE)
  mat <- xtv - mat
  if (order.g == 1) G <- 1/(1+exp(-speedadj*as.vector(mat)))
  if (order.g > 1) {
    prodmat <- mat[,1]
    for (i in 2:order.g) {
      prodmat <- prodmat*mat[,i]
    }
    G <- 1/(1+exp(-speedadj*prodmat))
  }
  if (as.zoo == TRUE) G <- zoo(G)
  if (verbose == TRUE) {
    if (as.zoo == TRUE) xtv <- zoo(xtv)
    return(list(G = G, xtv = xtv))
  }
  else return(G)
}	
