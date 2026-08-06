tvgarchObj <- function (par, fixed.par.g, y, order.g, xtv, opt, iter.fit.h, 
                        flag)
{
  n <- length(xtv)  
  s <- length(order.g)
  par.g <- c(fixed.par.g, par[1:(2*s+sum(order.g))])
  par.h <- par[-(1:(2*s+sum(order.g)))]
  G <- matrix(1, n, s)
  for (i in 1:s) {
    G[,i] <- 
      tv(speed = par.g[1+s+i], 
         location = par.g[(1+2*s+sum(order.g[1:i])-
                             order.g[i]+1):(1+2*s+sum(order.g[1:i]))], 
         xtv = xtv, opt = opt, order.g = order.g[i])
  }
  g <- par.g[1] + G %*% par.g[2:(s+1)]
  h <- garchxRecursion(pars = as.numeric(par.h), aux = iter.fit.h)
  if(flag != 2) {
    ll <- dnorm(x = y, mean = 0, sd = sqrt(g*h), log = TRUE)
  }
  if (flag == 0) {
    return(-ll)
  }
  if (flag == 1) {
    return(-sum(ll))
  }
}
