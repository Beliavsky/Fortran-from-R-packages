\name{riskSimul}
\alias{riskSimul}
\alias{riskSimul-package}
\docType{package}
\title{
Risk Quantification for Stock Portfolios under the T-Copula Model
}
\description{
This package can estimate the tail loss probabilities and conditional excess 
    for a stock portfolio. The log-returns are assumed to follow a t-copula model with generalized hyperbolic or t marginals. 
 
}
\details{
To simulate the tailloss probabilities of a portfolio for which the parameters of the t-copula model with generalized hyperbolic or t marginals
are available the following two functions can be used.

\code{\link{SISTCopula}()} is the name of the function that uses stratified importance sampling (SIS) to estimate 
a single or several tailloss probabilities and the corresponding conditional excess in a very efficient way.

\code{\link{NVTCopula}()} estimates the same quantities using naive simulation (without variance reduction).
}

\author{
  Wolfgang Hormann, Ismail Basoglu
}

\references{
I Basoglu, W Hormann. 2014. Efficient stratified sampling implementations in multiresponse 
simulation, in: Proceedings of the 2014 Winter Simulation Conference
A. Tolk, S. Y. Diallo, I. O. Ryzhov, L. Yilmaz, S. Buckley, and J. A. Miller, eds.

I Basoglu, W. Hormann, and H. Sak. 2013. Optimally Stratified Importance Sampling for Portfolio Risk
with Multiple Loss Thresholds. Optimization 62 (11): 1451-1471
} 



\examples{
R<- matrix(
c(1, 	0.554, 	0.632, 	0.419, 	0.400, 
  0.554,1, 		0.495, 	0.540, 	0.479,
  0.632,0.495, 	1, 		0.426, 	0.445,
  0.419,0.540, 	0.426, 	1, 		0.443,
  0.400,0.479, 	0.445, 	0.443, 	1),ncol=5)
  
pmg<- matrix(NA,ncol=5,nrow=5)  
colnames(pmg) <- c("lambda","alpha","beta","delta","mu")
pmg[1,] <- c(-0.602828, 8.52771, -0.533197, 0.014492, -0.000091)
pmg[2,] <- c(-1.331923, 2.72759, -2.573416, 0.019891, 0.001388)
pmg[3,] <- c(-1.602705, 3.26482, 1.456542, 0.035139, -0.001662)
pmg[4,] <- c(-1.131092, 15.13351, -1.722396, 0.014771, 0.001304)
pmg[5,] <- c(-0.955118, 31.14005, 0.896576, 0.015362, -0.000238)
 
portfo <- new.portfobj(nu=8.195,R=R,typemg="GH",parmg=pmg,c=rep(1,5),w=rep(0.2,5))

res1<- SISTCopula(n=10^4,npilot=c(10^3,3*10^3),portfobj=portfo,threshold=c(0.97,0.96,0.95,0.94),
                  stratasize=c(22,22),CEopt=FALSE,beta=0.75,mintype=0)
}










