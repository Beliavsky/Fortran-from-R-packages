\name{AsianCall}
\alias{AsianCall}
\title{
Calculates the Price, Delta and Gamma of an Asian Option
}
\description{
	Prices arithmetic average Asian Call options under geometric Brownian motion. It also estimates the sensitivities Delta and Gamma.
}
\usage{
AsianCall(T=1,d=12,K=100,r=0.05,sigma=0.2,S0=100,method=c("best","naive"),
          sampling=c("QMC","MC"),
          metpar=list(maxiter=100,tol=1.e-14,cvmethod="splitting"),
          sampar=list(nout=50,seq.type="korobov",n=2039,a=1487,
                      baker=TRUE,genmethod="pca"))
}
\arguments{
  \item{T}{
time to maturity (in years) 
}
  \item{d}{
  	number of control points
}
  \item{K}{
  	strike price
}
  \item{r}{
  	risk free interest rate
}
  \item{sigma}{
volatility
}
  \item{S0}{
starting stockprice 
}
\item{method}{selects the simulation method; \cr
method \code{"best"} uses a variance reduction method
based on effective control variates and conditional Monte Carlo and is very effective. \cr
method \code{"naive"} is mainly provided for comparison purposes.
}
\item{sampling}{
sampling \code{QMC} uses the Quasi Monte Carlo method Korobov lattice for the simulation.\cr
sampling \code{MC} uses standard Monte Carlo for the simulation.\cr 
}
\item{metpar}{list holding extra parameters related to the simulation method \cr
  For \code{method="best"}:\cr
   \code{maxiter=} maximal no of iterations for Newton method\cr
   \code{tol=} error tolerance for Newton method\cr
  for \code{sampling="QMC"}:\cr
      \code{cvmethod=c("splitting","direct")}  NOT necessary for method = "naive"\cr
  		  "splitting" ... estimates CV coefficients using lm with bootstrap\cr
         "direct"  ... estimates CV coefficients using lm and the full sample\cr
for \code{sampling="MC"}:\cr
    \code{np} ... sample size for pilot run for CV; NOT necessary for method = "naive"\cr
}

\item{sampar}{list holding several parameters related to the sampling method; \cr
         for \code{sampling="MC"} the list \code{sampar} only contains the total samplesize \code{n};\cr   \cr
		 for \code{sampling="QMC"} the list \code{sampar} contains the elements:
\describe{
\item{\code{nout}}{number of independent "randomized" copies of the Korobov lattice}
\item{\code{n}}{number of points of the Korobov lattice}
\item{\code{a}}{important constant for the construction of the Korobov lattice}

\item{\code{baker}}{\code{TRUE/FALSE}, indicates if Baker transform should be used for making the integrand periodic} 
\item{\code{genmethod}}{\code{= c("pca", "std","pcamain","lt","ltpca")},\cr note that for 
                \code{method=="naive"} only \code{genmethod=c("pca","std")} can be used.
                 \describe{
					\item{\code{genmethod="pca"}}{principal component analysis}
					\item{\code{genmethod="std"}}{standard}
					\item{\code{genmethod="pcamain"}}{use only first \code{dirnum} main directions of the PCA}
					\item{\code{genmethod="lt"}}{uses a transform for the first \code{dirnum}}
					\item{\code{genmethod="ltpca"}}{combination of \code{lt} with \code{pca}}
				}				 
			}
\item{\code{dirnum}}{number of main directions, only used for \code{genmethod="pcamain"} or \code{"lt"}}
}
}
}
\details{
Method \code{best} (see the reference Dingec and Hormann below)
is a very efficient simulation algorithm using multiple Control Variates and conditional MonteCarlo 
to calculate the the price, delta and gamma of Asian call options under geometric Brownian motion.
It is especially effective when \code{QMC} is selected as sampling method.
As QMC method Korobov Lattice rules are used. For good parameter values
see Table 1 of (L'Ecuyer, Lemieux).
}
\value{
	 returns a matrix holding the price and greeks.
	 The estimated Asian Call price and its estimated delta and gamma 
  form the first column vector, the respective 95 percent error bounds are given in the second column.
}
\references{
K. D. Dingec and W. Hormann. 
Improved Monte Carlo and Quasi-Monte Carlo Methods for the Price and the Greeks of Asian Options, 
Proceedings of the 2014 Winter Simulation Conference
A. Tolk, S. D. Diallo, I. O. Ryzhov, L. Yilmaz, S. Buckley, and J. A. Miller, eds.

L'Ecuyer, P., and C. Lemieux. 2000. Variance Reduction via Lattice Rules. Management Science 46 (9): 1214-1235.
	}
\author{
Kemal Dingec, Wolfgang Hormann
}
\seealso{
\code{\link{OptionPricing-package}}
}

\examples{
# standard settings for an efficient simulation using QMC and variance reduction
AsianCall(T=1,d=12,K=100,r=0.05,sigma=0.2,S0=100,method="best",
  sampling="QMC",metpar=list(maxiter=100,tol=1.e-14,cvmethod="splitting"),
  sampar=list(nout=50,n=2039,a=1487,baker=TRUE,genmethod="pca"))
					   
# efficient Monte Carlo version of the above simulation
AsianCall(T=1,d=12,K=100,r=0.05,sigma=0.2,S0=100,method="best",
  sampling="MC",metpar=list(maxiter=100,tol=1.e-14,np=1000),
  sampar=list(n=10^5))

# simple QMC version without variance reduction

AsianCall(T=1,d=12,K=100,r=0.05,sigma=0.2,S0=100,method="naive",
  sampling="QMC",
  sampar=list(nout=50,n=2039,a=1487,baker=TRUE,genmethod="pca"))
					  
# naive Monte Carlo version
AsianCall(T=1,d=12,K=100,r=0.05,sigma=0.2,S0=100,method="naive",
  sampling="MC",sampar=list(n=10^5))

					  

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
