\name{SISTCopula}
\alias{SISTCopula}
\alias{NVTCopula}
\alias{new.portfobj}
\title{ 
Efficient tail-loss probability and conditional excess estimation for t-copula model
}
 
\description{
Using stratified importance sampling (SIS) or naive simulation (NV)
the tail-loss probabilities and conditional excess values for several threshold values 
are estimated for a stock portfolio. 
The logreturns of the stocks are assumed to follow a t-copula model with generalized hyperbolic or t marginals. 

}

\usage{
SISTCopula(n=10^5,npilot=c(10^4,2*10^4),portfobj,threshold=c(0.95,0.9),
           stratasize=c(22,22),CEopt=FALSE,beta=0.75,mintype=-1)

NVTCopula(n=10^5,  portfobj, threshold=c(0.95,0.9))

new.portfobj(nu,R,typemg="GH",parmg,c=rep(1,dim(R)[1]),w=c/sum(c))
}  
\arguments{
  \item{n}{total sample size}
  \item{npilot}{size of one or several pilot runs, the sum of them should be smaller than \code{n}/2}
  \item{portfobj}{object of portfolio parameters}
  \item{threshold}{one or several threshold values (they should be ordered)}
  \item{stratasize}{a vector of length two holding the number of strata}
  \item{CEopt}{TRUE ... minimize the overall error of Conditional Exess estimates, otherwise of tail-loss estimates}
  \item{beta}{weight of maximal threshold value used for calculating the intermediate threshold used for selecting the IS density,
               only used when \code{length(threshold)>1}}
  \item{mintype}{only used when \code{length(threshold)>1}; 0 ... minimize mean square errors, -1 ... minimize relative MSE,
                -2 ... minimize the maximal error, -3 minimize the maximal relative error; a positive integer \code{j}
   				indicates that the variance of the estimate for the j-th threshold is minimized.}				
  \item{nu}{degrees of freedom of the t-copula}
  \item{R}{correlation matrix of the t-copula}
  \item{typemg}{type of the marginal distribution, \code{"GH"} generalized hyperbolic distribution, \code{"t"} t-distribution}
  \item{parmg}{matrix holding in its rows the parameters of the marginal distribution; 
              for the generalized hyperbolic distribution each row holds the parameters
			  lambda, alpha, beta, delta and mu;
			  for the t-distribution each row holds the parameters mu, sigma and nu (degrees of freedom).}
  \item{c}{scale factor vector of the portfolio}
  \item{w}{portfolio weights}
  
}
%\details{}  
\value{
For the case that the variable \code{threshold} contains only one value a matrix containing 
the results for the tail-loss probability in the first row and that of the conditional excess in the second row
is returned.

In the case that several threshold values are considered, a list consisting of the result matrices
for tail-loss probabilities and for conditional excess and the vector of the threshold values
is returned.
}

%\references{
%}

\author{
Ismail Basoglu,  Wolfgang Hormann
}

%\seealso{ }


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
 res1
 SISTCopula(n=10^4,npilot=c(10^3,3*10^3),portfobj=portfo,threshold=0.94,
            stratasize=c(22,22),CEopt=FALSE)

 NVTCopula(n=10^4,portfobj=portfo,threshold=c(0.97,0.96,0.95,0.94))
 NVTCopula(n=10^4,portfobj=portfo,threshold=0.94)

########
# example with t-marginals

R<- matrix(
c(1, 	0.551, 	0.636, 	0.421, 	0.398, 
  0.551,1, 	0.496, 	0.540, 	0.477,
  0.636,0.496, 	1, 	0.428, 	0.447,
  0.421,0.540, 	0.428, 	1, 	0.444,
  0.398,0.477, 	0.447, 	0.444, 	1),ncol=5)
   

pmg<- matrix(NA,ncol=3,nrow=5)  
colnames(pmg) <- c("mu","sigma","nu")
pmg[1,] <- c(-0.000258, 0.013769, 1.78)
pmg[2,] <- c(0.000794, 0.012166, 2.64)
pmg[3,] <- c(-0.000837, 0.019616, 3.25)
pmg[4,] <- c(0.001041, 0.009882, 2.67)
pmg[5,] <- c(-0.000104, 0.010812, 3.10)

portfo <- new.portfobj(nu=7.525,R=R,typemg="t",parmg=pmg,c=rep(1,5),w=rep(0.2,5))

res1<- SISTCopula(n=10^4,npilot=c(10^3,3*10^3),portfobj=portfo,threshold=c(0.97,0.96,0.95,0.94),
                  stratasize=c(22,22),CEopt=FALSE,beta=0.75,mintype=0)
res1
SISTCopula(n=10^4,npilot=c(10^3,3000),portfobj=portfo,threshold=0.94,stratasize=c(22,22))

NVTCopula(n=10^4,portfobj=portfo,threshold=c(0.97,0.96,0.95,0.94))
NVTCopula(n=10^4,portfobj=portfo,threshold=0.94)
}
\keyword{copula}