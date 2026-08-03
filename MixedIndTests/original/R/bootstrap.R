#'@title Function to perform multiplier bootstrap for tests of randomness or independence
#'
#'@description This function simulates Cramer-von Mises statistics using Gaussian multipliers.
#'
#'@param M          n x n x m vector  with MM[i,j] = 1(Xi <= Xj) and C=mean(M[,j]);
#'@param n          length of the series.
#'@param J          n x n vector for bootstrapping Sn.
#
#'
#'
#'@return \item{cvm_sim}{Simulated value of the Cramer-von Mises statistics}
#'@return \item{sn_sim}{simulated value of the Cramer-von Mises statistic Sn}
#'
#'@references B.R Nasri (2021). Tests of serial dependence for arbitrary distributions
#'
#'
#'@keywords internal
#'
#'@export
#'
#'
bootstrap = function(M,J,n){

cA = length(M)/n/n
xi = rnorm(n)

out0 = .C("statsim",
          as.double(M),
          as.double(J),
          as.double(xi),
          as.integer(n),
          as.integer(cA),
          stats=double(cA),
          sn = double(1),
          PACKAGE = "MixedIndTests"
          )

out = list(cvm=out0$stats,Sn = out0$sn)
return(out)
}
