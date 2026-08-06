tvObj <- function (par.g, fixed.par.g, xtv, opt, order.g, fixed.h, y, iter0, 
                   flag)
{
  n <- length(xtv)  
  s <- length(order.g)
  if (iter0 == FALSE) {
    par.g <- c(fixed.par.g, par.g)
  }
  G <- matrix(1, n, s)
  for (i in 1:s) {
    G[,i] <- 
      tv(speed = par.g[1+s+i], 
         location = par.g[(1+2*s+sum(order.g[1:i])-
                             order.g[i]+1):(1+2*s+sum(order.g[1:i]))], 
         xtv = xtv, opt = opt, order.g = order.g[i])
  }
  g <- par.g[1] + G %*% par.g[2:(s+1)]
  if(flag != 2) {
    phi <- y / sqrt(fixed.h)
    ll <- dnorm(x = phi, mean = 0, sd = sqrt(g), log = TRUE)
  }
  if (flag == 0) {
    return(-ll)
  }
  if (flag == 1) {
    return(-sum(ll))
  }
  if (flag == 2) {
    return(g)
  }
}
