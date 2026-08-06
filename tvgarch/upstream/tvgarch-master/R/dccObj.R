dccObj <- function (par.dcc, z, sigma2, flag)
{
  n <- nrow(z)
  m <- ncol(z)
  dcc <- matrix(0, n, m*(m-1)/2)
  if (flag != 2) ll <- numeric(n)
  R <- cor(z)
  Qbar <- (1 - par.dcc[1] - par.dcc[2]) * R
  Qt <- cor(z[1:100,])
  zt2 <- cor(z[1:100,])
  for (t in 1:n) {
    Qt <- Qbar + par.dcc[1] * zt2 + par.dcc[2] * Qt
    invQt <- diag( 1 / sqrt( diag(Qt) ) )
    Rt  <- invQt %*% Qt %*% invQt
    zt2 <- z[t,] %o% z[t,]
    if (flag != 2) ll[t] <- - 0.5 * (m * log(2*pi) + sum(log(sigma2[t,])) + 
                                       log(det(Rt)) + 
                                       sum((z[t,] %*% solve(Rt)) * z[t,])) 
    if (flag == 2) dcc[t,] <- Rt[lower.tri(Rt)]
  }
  
  if (flag == 0) return(-ll)
  if (flag == 1) return(-sum(ll))
  if (flag == 2) return(dcc)
}

