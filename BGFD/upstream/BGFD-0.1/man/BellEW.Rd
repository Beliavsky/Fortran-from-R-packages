\name{The Bell exponentiated Weibull distribution}
\alias{dBellEW}
\alias{pBellEW}
\alias{qBellEW}
\alias{rBellEW}
\alias{sBellEW}
\alias{hBellEW}
\alias{mBellEW}
\title{The Bell exponentiated Weibull distribution
}

\description{Density, probability, quantile function, random generation, survival function, hazard rate function and maximum likelihood estimates from the Bell exponentiated Weibull distribution.}

\usage{
dBellEW(x, alpha, beta, theta, lambda, log = FALSE)
pBellEW(x, alpha, beta, theta, lambda, log.p = FALSE, lower.tail = TRUE)
qBellEW(p, alpha, beta, theta, lambda, log.p = FALSE, lower.tail = TRUE)
rBellEW(n, alpha, beta, theta, lambda)
sBellEW(x, alpha, beta, theta, lambda, log.p = FALSE, lower.tail = TRUE)
hBellEW(x, alpha, beta, theta, lambda, log = FALSE,log.p = FALSE, lower.tail = TRUE)
mBellEW(x, alpha, beta, theta, lambda, method="B")
}

\arguments{
\item{x}{
A vector of (non-negative integer) quantiles.
}
\item{p}{
A vector of probablities.
}
\item{n}{
The number of random values to be generated under the Bell exponentiated Weibull distribution.
}
\item{lambda}{
The strictly positive parameter of the Bell distribution (\eqn{\lambda > 0}).
}
\item{alpha}{
The strictly positive scale parameter of the baseline exponentiated Weibull distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive shape parameter of the baseline exponentiated Weibull distribution (\eqn{\beta > 0}).
}
\item{theta}{
The strictly positive shape parameter of the baseline exponentiated Weibull distribution (\eqn{\theta > 0}).
}
\item{lower.tail}{if FALSE then 1-F(x) are returned and quantiles are computed 1-p.
}

\item{log}{if TRUE, probabilities p are given as log(p).
}
\item{log.p}{if TRUE, probabilities p are given for exp(p).
}
\item{method}{the procedure for optimizing the log-likelihood function after setting the intial values of the parameters and data values for which the Bell exponentiated Weibull distribution is fitted. It could be "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", or "SANN". "BFGS" is set as the default.
}
}

\details{The functions allow fitting the compounded Bell exponentiated Weibull distribution and evaluating the probability density function, cumulative distribution function, 
quantile function, random numbers, survival function, hazard rate function, and maximum likelihood estimates (MLEs) of the unknown parameters with standard error (SE) of 
the estimates. It also provides the goodness-of-fit measures such as the Akaike Information Criterion (AIC), Bayesian Information Criterion (BIC), 
the minimum value of the negative log-likelihood function, Anderson-Darling (A) test, Cramer-Von-Mises (W) test, Kolmogorov-Smirnov test, P-value and convergence status. 
}



\value{
dBellEW gives the (log) probability function.
pBellEW gives the (log) distribution function.
qBellEW gives the quantile function.
rBellEW generates random values.
sBellEW gives the survival function.
hBellEW gives the hazard rate function.
mBellEW gives the estimated parameters along with SE and goodness-of-fit measures.
}

\references{
Fayomi, A., Tahir, M. H., Algarni, A., Imran, M. and Jamal, F. (2022). A new useful exponential model
with applications to quality control and actuarial data.
Computational Intelligence and Neuroscience.
 
Nadarajah, S., Cordeiro, G. M. and Ortega, E. M. (2013). The exponentiated Weibull distribution: a survey. 
Statistical Papers, 54: 839-877.

Castellares, F., Ferrari, S. L., and Lemonte, A. J. (2018). On the Bell distribution and its
associated regression model for count data. Applied Mathematical Modelling, 56: 172-185.
}

\author{
Muhammad Imran and Michail Tsagris.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com} and Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{pBellEE}
}
}

\examples{
x<-rBellEW(20,2,1,1.2,0.2)
dBellEW(x,2,1,1.2,0.2)
pBellEW(x,2,1,1.2,0.2)
qBellEW(0.7,2,1,1.2,0.2)
sBellEW(x,2,1,1.2,0.2)
hBellEW(x,2,1,1.2,0.2)
mBellEW(x, 0.2,0.1,1.2,0.2, method="B")
}

