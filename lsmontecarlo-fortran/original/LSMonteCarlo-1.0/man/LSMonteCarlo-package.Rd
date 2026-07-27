\name{LSMonteCarlo-package}
\alias{LSMonteCarlo-package}
\alias{LSMonteCarlo}
\docType{package}
\title{
American options pricing with Least Squares Monte Carlo method
}
\description{
The package compiles functions that calculate prices of American put options with Least Squares Monte Carlo method. The option types are plain vanilla American put, Asian American put, and Quanto American put. The pricing algorithms include variance reduction techniques such as Antithetic Variates and Control Variates. Additional functions are given to derive "price surfaces" at different volatilities and strikes, create 3-D plots, quickly generate Geometric Brownian motion, and calculate prices of European options with Black & Scholes analytical solution.
}
\details{
\tabular{ll}{
Package: \tab LSMonteCarlo\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-09-20\cr
License: \tab GPL 3\cr
}
The Least Squares Monte Carlo is an approach developed to approximate the value of American options. It combines regression modeling and Monte Carlo simulation. The key feature of this method is estimation of the conditional expectation of the future pay-offs by a regression model (for details see Longstaff & Schwartz, 2000). The main pricing functions employing this method in the package are: \code{AmerPutLSM}, \code{AsianAmerPutLSM}, and \code{QuantoAmerPutLSM}. Pricing functions that include variance reduction methods are: \code{AmerPutLSM_AV}, \code{QuantoAmerPutLSM_AV} (Antithetic Variates) and \code{AmerPutLSM_CV} (Control Variates, with Black & Scholes solution for European put used as the control). All these functions are based on Geometric Brownian motion as a price process. They can be used with tailored \code{summary}, \code{print}, and \code{price} functions. The "price surfaces" at different volatilities and strikes can be derived using the functions \code{AmerPutLSMPriceSurf}, \code{AsianAmerPutLSMPriceSurf}, and \code{QuantoAmerPutLSMPriceSurf}, and plotted with tailored \code{plot} function. For general reading on option pricing with Monte Carlo methods see Glasserman (2004).
}
\author{
Mikhail A. Beketov

Maintainer: Mikhail A. Beketov <mikhail.beketov@gmx.de>
}
\references{
Glasserman, P. 2004. Monte Carlo Methods in Financial Engineering. Springer.

Longstaff, F.A., and E.S. Schwartz. 2000. Valuing american option by simulation: A simple least-squared approach. The Review of Financial Studies. 14:113-147.

}

\keyword{ 
Quantitative Finance
Option pricing
Monte Carlo
Regression 
}
\seealso{
Functions: \code{\link{AmerPutLSM}}, 
\code{\link{AmerPutLSM_AV}}, 
\code{\link{AmerPutLSM_CV}}, 
\code{\link{AsianAmerPutLSM}}, 
\code{\link{QuantoAmerPutLSM}}, and 
\code{\link{QuantoAmerPutLSM_AV}}.
}
\examples{
Put<-AmerPutLSM(Spot=14.2, Strike=16.5, n=200, m=50)
summary(Put)
price(Put)
plot(AmerPutLSMPriceSurf(vols = (seq(0.1, 1.5, 0.2)), n=200, m=10, 
strikes = (seq(0.5, 1.9, 0.2))), color = divPalette(150, "RdBu"))
}
