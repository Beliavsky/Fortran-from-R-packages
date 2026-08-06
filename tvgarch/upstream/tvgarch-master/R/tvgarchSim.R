tvgarchSim <- function(n, order.g = 1, order.h = c(1,1,0), intercept.g = 1.2, 
                       size = 5, speed = 25, location = 0.5, xtv = NULL, 
                       intercept.h = 0.2, arch = 0.1, garch = 0.8, asym = NULL, 
                       xreg = NULL, opt = 0, as.zoo = TRUE, verbose = FALSE, 
                       innovations = NULL)
{
  m <- 1
  names.ID <- "y"
  if (!is.null(order.g)) {
    order.g <- as.matrix(order.g)
    s <- length(order.g)
    if (is.null(intercept.g) || is.null(size) || is.null(speed) || is.null(location)) stop("At least one value of the parameters in the TV component is missing.")
    size <- as.matrix(size)
    if (s != length(size)) stop("Mismatch between the number of transition functions, s, and the number of size parameters.")
    speed <- as.matrix(speed)
    if (s != length(speed)) stop("Mismatch between the number of transition functions, s, and the number of speed parameters.")
    location <- as.matrix(location)
    if (sum(order.g) != length(location)) stop("Mismatch between order.g and the number of location parameters.")
  }
  '
    Constructing the g component and standardising residuals by sqrt(g)
  '
  if (!is.null(order.g)) { 
    if (is.null(xtv)) xtv <- (1:n)/n
    G <- matrix(1, n, s)
    for (j in 1:s) {
      G[,j] <- 
        tv(speed = speed[j], 
           location = location[(sum(order.g[1:j])-
                                  order.g[j]+1):sum(order.g[1:j])], xtv = xtv, 
           opt = opt, order.g = order.g[j])
    }
    g <- intercept.g + G %*% size
    colnames(g) <- "g"
  }
  if (is.null(order.g)) g <- matrix(1, n, m)
  '
    Constructing the h component
  '
  hSim <- garchxSim(n = n, intercept = intercept.h, arch = arch, garch = garch, 
                    asym = asym, xreg = xreg, innovations = innovations, 
                    as.zoo = FALSE, verbose = TRUE)
  h <- as.matrix(hSim[,"sigma2"])
  colnames(h) <- "h"
  innovations <- as.matrix(hSim[,"innovations"])
  colnames(innovations) <- "innovations"
  '
    Output
  '
  sigma2 <- h*g
  colnames(sigma2) <- "sigma2"
  y <- sqrt(sigma2)*innovations
  colnames(y) <- names.ID
  Ineg <- as.matrix(as.numeric(innovations < 0))
  colnames(Ineg) <- "Ineg"
  if (verbose == TRUE) {
    result <- cbind(y, sigma2, g, h, innovations, Ineg)
  }
  if(verbose == FALSE){
    result <- y
  }
  result <- as.matrix(result)
  if (as.zoo == TRUE) {
    result <- as.zoo(result)
  }
  return(result)
}
