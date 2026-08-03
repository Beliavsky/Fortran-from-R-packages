#'@title Simulation of a copula-based time series
#'
#'@description This function simulates a Markovian time series (p-Markov for the Farlie-Gumbel-Morgenstern copula) with uniform margins using a copula family for the joint distribution of U[t], U[t-1].
#'
#'@param family  "ind", "tent", "gaussian", "t" , "clayton", "fgm", "frank", "gumbel", joe" , "plackett"
#'@param n       length of the time series
#'@param tau     Kendall's tau of the copula family
#'@param param   extra copula parameter: for "fgm",  param is the dimension of the copula; for "t", param = nu
#'
#'@return \item{U}{Simulated time series}
#'
#'@author Bouchra R. Nasri January 2021
#'@references B.R Nasri (2022). Tests of serial dependence for arbitrary distributions
#'
#'@examples
#' U = SimCopulaSeries("fgm",100,0.2, 3) # for the FGM, |tau|<= 2/9
#'@export


SimCopulaSeries=  function(family,n,tau=0,param=NULL)
{

  switch(family,
          "tent"={
                  U  = numeric(n)
                  U[1] = runif(1)
                  for(i in 2:n){ U[i] = 1.999*min(U[i-1],1-U[i-1]) }
                 },

    "clayton" = {
                  if(tau==0){U=runif(n)}
      else{
                   theta = copula::iTau(copula::claytonCopula(),  tau)
                   f = function(u,w){
                           p = ( 1+ u^(-theta) *(-1 + w^(-theta/(1+theta))))^(-1/theta)
                         return(p)
                       }
                   U=numeric(n)
                   U[1] = runif(1);
                   for(i in 2:n)
                      {
                        U[i] = f(U[i-1],runif(1))
                   }
      }
               },

     "gaussian"={
        rho = copula::iTau(copula::normalCopula(),  tau)
        cte =sqrt(1-rho^2)
        Z=numeric(n)
        Z[1] = rnorm(1)
         for(i in 2:n){
            Z[i] = rho*Z[i-1]+ cte*rnorm(1)
         }
        U = pnorm(Z)
},
    "t" = {
            nu = param$nu;
            rho = copula::iTau(copula::normalCopula(),  tau)
            U = numeric(n)
            U[1] = runif(1)
            x = qt(U[1],nu)

           Omega = sqrt(1-rho^2)
           sn = 1/sqrt(nu+1);
           for( i in 2:n){
               mu = rho*x
               v = runif(1)
               y = qt(v,nu+1)
               U[i] = pt( mu+y*Omega*sqrt(nu+x^2)*sn,nu)
               x = qt(U[i],nu)
                 }
            },

    "frank"={
      if(tau==0){U=runif(n)}
      else{
        theta =  copula::iTau(copula::frankCopula(),  tau);
           f = function(u,w) { p = log( 1 -w*(1-exp(-theta))/(exp(-theta*u)+w*(1-exp(-theta*u))) )* (-1/theta)
           return(p)}

           U=numeric(n)
           U[1] = runif(1);
           for(i in 2:n){  U[i] = f(U[i-1],runif(1))}
      }
        },

   "fgm"={
     p = param
     theta = copula::iTau(copula::fgmCopula(),tau)
        U =numeric(n)
         U[1:p-1]=runif(p-1);
          for(j in p:n)
          {
            w = runif(1);
             A = theta*prod(1-2*U[(j-p+1):(j-1)]);
            U[j] = 2*w/(1+A+sqrt((1+A)^2-4*A*w));
          }
        },
    "ind" ={ U=runif(n) },

"plackett" = {
  if(tau==0){U=runif(n)}
  else{
               theta = copula::iTau(copula::plackettCopula(),  tau)
               U = numeric(n)
               U[1] = runif(1)
               f = function(u,w){
                                  v = w*(1-w)
                                  A = theta+(theta-1)^2*v
                                  C = v* (1+(theta-1)*u)^2
                                  B = 2*(theta-1)*v* (1-(theta+1)*u)-theta
                                  p = 0.5*(-B + sign(w-0.5)*sqrt(B^2-4*A*C))/A
                          return(p)
               }
               for(i in 2:n){  U[i] = f(U[i-1],runif(1))}
  }
},
    "gumbel" ={
      theta = copula::iTau(copula::gumbelCopula(),  tau)
      cop <- copula::gumbelCopula(theta, dim = 2)
      U = numeric(n)
      U[1] = runif(1)
      f = function(u,w){ return(copula::cCopula(c(u,w), cop, indices = 1:2, inverse = TRUE))      }
      for(i in 2:n){  U[i] = f(U[i-1],runif(1))}

    },
"joe" ={
  theta = copula::iTau(copula::joeCopula(),  tau)
  cop <- copula::joeCopula(theta, dim = 2)
  U = numeric(n)
  U[1] = runif(1)
  f = function(u,w){ return(copula::cCopula(c(u,w), cop, indices = 1:2, inverse = TRUE))      }
  for(i in 2:n){  U[i] = f(U[i-1],runif(1))}

}
     )

return(U)
}

