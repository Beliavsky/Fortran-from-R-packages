\name{The Bell Fisk distribution}
\alias{dBellF}
\alias{pBellF}
\alias{qBellF}
\alias{rBellF}
\alias{sBellF}
\alias{hBellF}
\alias{mBellF}
\title{The Bell Fisk distribution
}

\description{Density, probability, quantile function, random generation, survival function, hazard rate function and maximum likelihood estimates from the Bell Fisk distribution.}


\usage{
dBellF(x, a, b, lambda, log = FALSE)
pBellF(x, a, b, lambda, log.p = FALSE, lower.tail = TRUE)
qBellF(p, a, b, lambda, log.p = FALSE, lower.tail = TRUE)
rBellF(n, a, b, lambda)
sBellF(x, a, b, lambda, log.p = FALSE, lower.tail = TRUE)
hBellF(x, a, b, lambda, log = FALSE,log.p = FALSE, lower.tail = TRUE)
mBellF(x, a, b, lambda, method="B")
}

\arguments{
\item{x}{
A vector of (non-negative integer) quantiles.
}
\item{p}{
A vector of probablities.
}
\item{n}{
The number of random values to be generated under the Bell Fisk distribution.
}
\item{lambda}{
The strictly positive parameter of the Bell distribution (\eqn{\lambda > 0}).
}
\item{a}{
The strictly positive scale parameter of the baseline Fisk distribution (\eqn{a > 0}).
}
\item{b}{
The strictly positive shape parameter of the baseline Fisk distribution (\eqn{b > 0}).
}

\item{lower.tail}{if FALSE then 1-F(x) are returned and quantiles are computed 1-p.
}

\item{log}{if TRUE, probabilities p are given as log(p).
}
\item{log.p}{if TRUE, probabilities p are given for exp(p).
}
\item{method}{the procedure for optimizing the log-likelihood function after setting the intial values of the parameters and data values for which the Bell Fisk distribution is fitted. It could be "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", or "SANN". "BFGS" is set as the default.
}
}

\details{The functions allow fitting the compounded Bell Fisk distribution and evaluating the probability density function, cumulative distribution function, 
quantile function, random numbers, survival function, hazard rate function, and maximum likelihood estimates (MLEs) of the unknown parameters with standard error (SE) of 
the estimates. It also provides the goodness-of-fit measures such as the Akaike Information Criterion (AIC), Bayesian Information Criterion (BIC), 
the minimum value of the negative log-likelihood function, Anderson-Darling (A) test, Cramer-Von-Mises (W) test, Kolmogorov-Smirnov test, P-value and convergence status. 
}


\value{
dBellF gives the (log) probability function.
pBellF gives the (log) distribution function.
qBellF gives the quantile function.
rBellF generates random values.
sBellF gives the survival function.
hBellF gives the hazard rate function.
mBellF gives the estimated parameters along with SE and goodness-of-fit measures.
}

\references{
Fayomi, A., Tahir, M. H., Algarni, A., Imran, M. and Jamal, F. (2022). A new useful exponential model
with applications to quality control and actuarial data. Computational Intelligence and Neuroscience, 2022.

Kleiber, C. and Kotz, S. (2003). Statistical size distributions in economics and actuarial sciences. 
John Wiley & Sons.

Castellares, F., Ferrari, S. L., and Lemonte, A. J. (2018). On the Bell distribution and its
associated regression model for count data. Applied Mathematical Modelling, 56:172-185
}

\author{
Muhammad Imran and Michail Tsagris.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com} and Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{pBellW}
}
}

\examples{
x<-rBellF(20,2,1,1.7)
x
dBellF(x,2,1,1)
pBellF(x,2,1,1)
qBellF(0.7,2,1,1)
sBellF(x,2,1,1.2)
hBellF(x,2,1,1.2)
mBellF(x, 0.2,1.1,0.7, method="B")
}

