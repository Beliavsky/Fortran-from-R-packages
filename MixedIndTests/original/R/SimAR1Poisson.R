#'@title Simulation of a AR(1) Poisson process
#'
#'@description Conditionally on the past, X[t] is Poisson with lambda[t] = a+bX[t-1]
#'
#'@param param     Param[1] = a>0, param[2] = b, 0<=b <1 (for stationarity)
#'@param n         Length of the series.
#
#'
#'
#'@return \item{X}{Simulated series}
#'
#'
#'@examples
#' data <- SimAR1Poisson(c(5,0.4),500)
#'@export
#'
#'


SimAR1Poisson = function(param,n)
{
n0 = n+100
X0 = 0*c(1:n0);



X0[1] = rpois(1,param[1]);
for( i in c(2:n0)){
  lambda = param[1]+param[2]*X0[i-1]
X0[i] = rpois(1,lambda);
}
X = X0[101:n0];
return(X)
}
