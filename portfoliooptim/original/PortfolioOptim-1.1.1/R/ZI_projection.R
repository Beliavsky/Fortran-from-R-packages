#'Computes \cr
#'the solution of the linear program
#'\eqn{\min c^{T} x}
#'  for \eqn{Ax = b}, \eqn{x \ge 0}
#'
#'   which fulfills the condition \eqn{|B(x - xhat)|^2 \to \min}
#'
#'@param clin g
#'@param Amat g
#'@param bmat g
#'@param xhat g
#'@param B g 
#'@param maxiter g
#'@param tol g
#'@param k g
#'
#'@return list
#'@keywords internal 

.ZI_projection <- function (clin, Amat, bmat, xhat, B, maxiter, tol, k)
{
  # fixed values of some parameters  
  tol1 = 1e-3
  delta = 0.5
  p_pow = 0.8
  tol_p = sqrt(tol)
  
  
  n = ncol(Amat)
  m = nrow(Amat)
  bmat = cbind(bmat)
  xhat = cbind(xhat)
  clin = cbind(clin)
  en = cbind(rep(1, n))
  em = cbind(rep(1, m))
  
  x_old = en + xhat
  x_old[x_old < 1] <- 1
  y_old = em
   
  rho_old = pmax(1, max(rbind(abs(t(Amat)%*%y_old - clin), abs( - Amat %*% x_old + bmat))) + 1) + delta
  s_old = t(B) %*% (B %*% (x_old-xhat))
  
  kap = min(1, 1/max(x_old * s_old), 1/max(x_old))  
  
  s_old = rho_old^p_pow *  kap * s_old
  s_old[s_old < tol1] <- kap
  z_old = rho_old^p_pow * kap * em
  
  F_old = .F_func(rho_old, x_old, y_old, s_old, z_old, Amat, bmat, clin, p_pow, B, xhat, kap )
  
  eta = max(abs(F_old))/rho_old
  bet = (eta + 1)/2
  ind =0
  
  while (rho_old > tol)
  {
    ind = ind+1
    
    x_old = as.numeric(x_old)  
    y_old = as.numeric(y_old)  
    z_old = as.numeric(z_old)  
    s_old = as.numeric(s_old)  
    
    #if ((ind %% 20) == 0) print(rho_old)
    if (max(abs(F_old)) == 0)
    {
      x_new = x_old
      y_new = y_old
      s_new = s_old
      z_new = z_old
    }
    else
    {
      x_diag = .make_diag(x_old)
      y_diag = .make_diag(y_old)
      rho_pow = rho_old^p_pow
      if (n > m) {
        G_inv = solve(.make_diag(s_old) + rho_pow * kap * x_diag %*% (t(B) %*% B))
        Gk =  ((y_diag) %*% Amat) %*% G_inv
        Hk = .make_diag(z_old) + rho_pow * kap * y_diag + Gk %*% (x_diag %*% t(Amat))
        Vl = (rho_old*em - y_diag %*% (rho_pow * kap * y_old + Amat %*% x_old - bmat)
              - Gk %*% (rho_old*en - x_diag %*%(rho_pow * kap* t(B) %*% (B %*% (x_old -xhat)) - t(Amat)%*% y_old +clin)))
        dy = solve(Hk, Vl)
        dx = G_inv %*% ((x_diag %*% t(Amat)) %*% dy + rho_old *en
                        - x_diag %*% (rho_pow * kap * t(B) %*% (B %*% (x_old - xhat)) - t(Amat) %*% y_old +clin))
      }
      else {
        N_inv = solve(.make_diag(z_old) + rho_pow * kap * y_diag)
        Nk = (x_diag %*% t(Amat)) %*% N_inv
        Mk = .make_diag(s_old) + rho_pow * kap* t(B) %*% (B %*% x_diag) + Nk %*% (y_diag %*% Amat)
        Wl = (rho_old*en - x_diag %*% (rho_pow * kap *t(B) %*% (B %*% (x_old - xhat)) - t(Amat) %*% y_old + clin)
              + Nk %*% (rho_old * em - y_diag %*% (Amat %*% x_old + rho_pow * kap * y_old -bmat)))
        dx = solve(Mk, Wl)
        dy = N_inv %*% 	(-(y_diag %*% Amat) %*% dx + rho_old * em
                         - y_diag %*% (Amat %*% x_old + rho_pow * kap * y_old -bmat))
      }
      ds = rho_pow * kap * t(B) %*%(B %*% dx) - t(Amat) %*% dy + rho_pow * kap * t(B) %*% (B %*% (x_old - xhat)) - t(Amat) %*% y_old - s_old +clin
      dz = Amat %*% dx + rho_pow * kap *dy + Amat %*% x_old + rho_pow * kap * y_old - z_old - bmat
      xt <- x_old/dx
      if (min(xt) >=0) l1 <- Inf else l1 <- min(-xt[xt<0])
      yt <- y_old/dy
      if (min(yt) >=0) l2 <- Inf else  l2 <- min(-yt[yt<0])
      st <- s_old/ds
      if (min(st) >=0) l3 <- Inf else  l3 <- min(-st[st<0])
      zt <- z_old/dz
      if (min(zt) >=0) l4 <- Inf else l4 <- min(-zt[zt<0])
      
      alpha = min(1,l1,l3,l4)
      
      # "Armijo rule" for lambda
      # set parameters  sigma = 0.001, b1 = 0.9
      # the rule is lambda_j  = alpha * b1^j   where
      # |F(x+lambda_j *dx)| <= (1 -sigma * lambda_j) |F(x)|
      sigma = 0.0001
      b1 = 0.9
      
      lkn = alpha
      Fdo = .F_func(rho_old, x_old, y_old, s_old, z_old,  Amat, bmat, clin, p_pow, B, xhat, kap )
      Fdo_norm = max(abs(Fdo))
      Fdn = .F_func(rho_old, x_old + lkn *dx, y_old + lkn *dy, s_old + lkn *ds, z_old + lkn *dz,  Amat, bmat, clin, p_pow, B, xhat, kap)
            
      while (max(abs(Fdn)) > (1-sigma*lkn) * Fdo_norm)
      {
        lkn = lkn * b1
        Fdn = .F_func(rho_old, x_old + lkn *dx, y_old + lkn *dy, s_old + lkn *ds, z_old + lkn *dz,  Amat, bmat, clin, p_pow, B, xhat, kap)
        
      }
      
      x_new = x_old + lkn * dx
      y_new = y_old + lkn * dy
      s_new = s_old + lkn * ds
      z_new = z_old + lkn * dz
    }
    # "Armijo rule" for gamma
    # set parameter    b2 = 0.9
    # the rule is gamma_j  =  b1^j   where
    # |F((1-gamma_j)*rho_old, x_new)| <= bet * (1 -gamma_j)* rho_old
    b2 = 0.9
    
    gkn = b2
    rho_new = (1-gkn)*rho_old
    ll = 0
    Fgn = .F_func(rho_new, x_new, y_new, s_new, z_new,  Amat, bmat, clin, p_pow, B, xhat, kap)
     
    while (max(abs(Fgn)) > bet *rho_new)
    {
      ll = ll+1
      gkn = gkn *b2
      rho_new = (1-gkn)*rho_old
      Fgn = .F_func(rho_new, x_new, y_new, s_new, z_new,  Amat, bmat, clin, p_pow, B, xhat, kap)
    }
    
    x_old = x_new
    y_old = y_new
    s_old = s_new
    z_old = z_new
    rho_old = rho_new
    
    if (ind > maxiter){
			cat(c("Maximal number of iteration has been exceeded. \n"))
			break
	}
	x_port = x_old[1:k]
	if (length(x_port[x_port < tol_p]) >= (k-1)){
			cat(c("Solution is a maximal corner solution, i.e. whole investment goes into one asset. \n"))
			break
	} 
  }
   
  return(list(xproj = x_old, yproj = y_old,  fmin = max(abs(Fgn)) ))
}
