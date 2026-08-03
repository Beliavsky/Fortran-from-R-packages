#'@title Kendall's tau and Spearman's rho statistics for testing randomness in a univariate time series
#'
#'@description This function computes  Kendall's tau and Spearman's rho statistics for tests of randomness in a time series with arbitrary distribution for pairs (X[i],X[i+k]), k=1:lags
#'
#'@param x      Time series
#'@param lag    Number of lags
#'@param graph  Set to TRUE for a dependogram for Kendall's tau and Spearman;s rho
#'
#'
#'@return \item{stat}{List of Kendall's tau and Spearman's rho statistics from multilinear copula, and test combinations LB}
#'@return \item{pvalue}{ P-values for the tests statistics}
#'
#'@references B.R Nasri (2022). Tests of serial dependence for arbitrary distributions
#'
#'@examples
#' out <-EstDepSerial(SimAR1Poisson(c(5,0.4),100),10)






EstDepSerial = function(x,lag,graph=FALSE)
{



n = length(x);
tau = 0*c(1:lag)
rho = 0*c(1:lag)
for(j in 1:lag)
{

  out0 = stat_dep_ser(x,j);
  tau[j] = out0$tau;
  rho[j] = out0$rho;

}
varH = out0$s2;
zrho= sqrt(n)*(rho);
pvalrho = 200*pnorm(-abs(zrho));
LBrho = sum(zrho^2);
pvalLBrho = 100*(1-pchisq(LBrho,lag));


z = sqrt(n)*tau;

ztau =  0.5*z/varH;
pvaltau = 200*pnorm(-abs(ztau));
LBtau = sum(ztau^2);
pvalLBtau = 100*(1-pchisq(LBtau,lag));

out0 = list(tau=tau,rho=rho,varH=varH,
           LBtau=LBtau,LBrho=LBrho,pvaltau=pvaltau,
           pvalrho=pvalrho,pvalLBtau=pvalLBtau,pvalLBrho=pvalLBrho,
           subsets=c(1:lag))
if(graph){AutoDep(out0)}
return(out0)
}


