\name{The Bell Weibull distribution}
\alias{dBellW}
\alias{pBellW}
\alias{qBellW}
\alias{rBellW}
\alias{sBellW}
\alias{hBellW}
\alias{mBellW}
\title{The Bell Weibull distribution
}

\description{Density, probability, quantile function, random generation, survival function, hazard rate function and maximum likelihood estimates from the Bell Weibull distribution.}
\usage{
dBellW(x, alpha, beta, lambda, log = FALSE)
pBellW(x, alpha, beta, lambda, log.p = FALSE, lower.tail = TRUE)
qBellW(p, alpha, beta, lambda, log.p = FALSE, lower.tail = TRUE)
rBellW(n, alpha, beta, lambda)
sBellW(x, alpha, beta, lambda, log.p = FALSE, lower.tail = TRUE)
hBellW(x, alpha, beta, lambda, log = FALSE,log.p = FALSE, lower.tail = TRUE)
mBellW(x, alpha, beta, lambda, method="B")
}

\arguments{
\item{x}{
A vector of (non-negative integer) quantiles.
}
\item{p}{
A vector of probablities.
}
\item{n}{
The number of random values to be generated under the Bell Weibull distribution.
}
\item{lambda}{
The strictly positive parameter of the Bell distribution (\eqn{\lambda > 0}).
}
\item{alpha}{
The strictly positive scale parameter of the baseline Weibull distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive shape parameter of the baseline Weibull distribution (\eqn{\beta > 0}).
}
\item{lower.tail}{if FALSE then 1-F(x) are returned and quantiles are computed 1-p.
}

\item{log}{if TRUE, probabilities p are given as log(p).
}
\item{log.p}{if TRUE, probabilities p are given for exp(p).
}
\item{method}{the procedure for optimizing the log-likelihood function after setting the intial values of the parameters and data values for which the Bell Weibull distribution is fitted. It could be "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", or "SANN". "BFGS" is set as the default.
}
}

\details{The functions allow fitting the compounded Bell Weibull distribution and evaluating the probability density function, cumulative distribution function, 
quantile function, random numbers, survival function, hazard rate function, and maximum likelihood estimates (MLEs) of the unknown parameters with standard error (SE) of 
the estimates. It also provides the goodness-of-fit measures such as the Akaike Information Criterion (AIC), Bayesian Information Criterion (BIC), 
the minimum value of the negative log-likelihood function, Anderson-Darling (A) test, Cramer-Von-Mises (W) test, Kolmogorov-Smirnov test, P-value and convergence status. 
}

\value{
dBellW gives the (log) probability function.
pBellW gives the (log) distribution function.
qBellW gives the quantile function.
rBellW generates random values.
sBellW gives the survival function.
hBellW gives the hazard rate function.
mBellW gives the estimated parameters along with SE and goodness-of-fit measures.
}

\references{

Fayomi, A., Tahir, M. H., Algarni, A., Imran, M. and Jamal, F. (2022). A new useful exponential model
with applications to quality control and actuarial data. Computational Intelligence and Neuroscience.

Castellares, F., Ferrari, S. L., and Lemonte, A. J. (2018). On the bell distribution and its
associated regression model for count data. Applied Mathematical Modelling, 56:172-185.
}

\author{
Muhammad Imran and Michail Tsagris.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com} and Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{pBellE}
}
}

\examples{
x<-rBellW(20,2,1,0.7)
dBellW(x,2,1,1)
pBellW(x,2,1,1)
qBellW(0.7,2,1,1)
sBellW(x,2,1,1.2)
hBellW(x,2,1,1.2)
mBellW(x, 0.2,0.1,1.2, method="B")
}

