#'Auxiliary function used by .ZI_projection
#'
#'@param rho g
#'@param x g
#'@param y g
#'@param s g
#'@param z g
#'@param A g
#'@param b g
#'@param cl g
#'@param p g
#'@param B g
#'@param xe g
#'@param ka g
#'@keywords internal 


.F_func <- function(rho, x, y, s, z, A, b, cl, p, B, xe, ka)
{
  n = ncol(A)
  m = nrow(A)
  x = cbind(x)
  y = cbind(y)
  s = cbind(s)
  z = cbind(z)
  b = cbind(b)
  cl = cbind(cl)
  xe = cbind(xe)


  F1 = x*s -rho*cbind(rep(1, n))
  F2 = y*z -rho*cbind(rep(1, m))
  F3 = s + t(A) %*% y - cl -rho^p * ka * t(B) %*% (B %*% (x -xe))
  F4 = z - A %*% x +b - rho^p* ka* y

  FF = rbind(F1, F2, F3, F4)

  return(FF)

}
