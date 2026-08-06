garchObj <- function (par.h, xreg, order.h, fixed.g, y, flag)
{
  n <- length(y)  
  if (order.h[1] != 0) par.garch <- par.h[(order.h[2]+2):(sum(order.h[1:2])+1)]
  if (order.h[1] == 0) par.garch <- NULL
  if (order.h[3] != 0) par.asym <- par.h[(sum(order.h[1:2])+2):(sum(order.h)+1)]
  if (order.h[3] == 0) par.asym <- NULL
  if (!is.null(xreg)) par.xreg <- xreg %*% par.h[(sum(order.h)+2):length(par.h)]
  if (is.null(xreg)) par.xreg <- NULL
  phi <- y / sqrt(fixed.g)
  h <- garchxSim(n = n, intercept = par.h[1], arch = par.h[2:(order.h[2]+1)], 
                 garch = par.garch, asym = par.asym, xreg = par.xreg, 
                 innovations = phi, as.zoo = FALSE, verbose = TRUE)[,"sigma2"]
  if (flag != 2) ll <- dnorm(x = phi, mean = 0, sd = sqrt(h), log = TRUE)    
  if (flag == 0) return(-ll)
  if (flag == 1) return(-sum(ll))
  if (flag == 2) return(h)
}
