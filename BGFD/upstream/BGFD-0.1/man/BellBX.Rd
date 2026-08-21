\name{The Bell Burr-X distribution}
\alias{dBellBX}
\alias{pBellBX}
\alias{qBellBX}
\alias{rBellBX}
\alias{sBellBX}
\alias{hBellBX}
\alias{mBellBX}
\title{The Bell Burr-X distribution}

\description{Density, probability, quantile function, random generation, survival function, hazard rate function and maximum likelihood estimates from the Bell Burr-X distribution.}

\usage{
dBellBX(x, a, lambda, log = FALSE)
pBellBX(x, a, lambda, log.p = FALSE, lower.tail = TRUE)
qBellBX(p, a, lambda, log.p = FALSE, lower.tail = TRUE)
rBellBX(n, a, lambda)
sBellBX(x, a, lambda, log.p = FALSE, lower.tail = TRUE)
hBellBX(x, a, lambda, log = FALSE, log.p = FALSE, lower.tail = TRUE)
mBellBX(x, a, lambda, method="B")
}

\arguments{
\item{x}{
A vector of (non-negative integer) quantiles.
}
\item{p}{
A vector of probablities.
}
\item{n}{
The number of random values to be generated under the Bell Burr-X distribution.
}
\item{lambda}{
The strictly positive parameter of the Bell distribution (\eqn{\lambda > 0}).
}
\item{a}{
The strictly positive shape parameter of the baseline Burr-X distribution (\eqn{a > 0}).
}
\item{lower.tail}{if FALSE then 1-F(x) are returned and quantiles are computed 1-p.
}

\item{log}{if TRUE, probabilities p are given as log(p).
}
\item{log.p}{if TRUE, probabilities p are given for exp(p).
}
\item{method}{the procedure for optimizing the log likelihood function after setting the intial values of the parameters and data values for which the BellBurr-X distribution is fitted. It could be "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", or "SANN". "BFGS" is set as the default.
}
}

\details{The functions allow fitting the compounded Bell-Burr-X distribution and evaluating the probability density function, cumulative distribution function, 
quantile function, random numbers, survival function, hazard rate function, and maximum likelihood estimates (MLEs) of the unknown parameters with standard error (SE) of 
the estimates. It also provides the goodness-of-fit measures such as the Akaike Information Criterion (AIC), Bayesian Information Criterion (BIC), 
the minimum value of the negative log-likelihood function, Anderson-Darling (A) test, Cramer-Von-Mises (W) test, Kolmogorov-Smirnov test, P-value and convergence status. 
}

\value{
dBellBX gives the (log) probability function.
pBellBX gives the (log) distribution function.
qBellBX gives the quantile function.
rBellBX generates random values.
sBellBX gives the survival function.
hBellBX gives the hazard rate function.
mBellBX gives the estimated parameters along with SE and goodness-of-fit measures.
}

\references{
Fayomi, A., Tahir, M. H., Algarni, A., Imran, M. and Jamal, F. (2022). A new useful exponential model
with applications to quality control and actuarial data.
Computational Intelligence and Neuroscience.

Kleiber, C. and Kotz, S. (2003). Statistical size distributions in economics and actuarial sciences. 
John Wiley & Sons.

Castellares, F., Ferrari, S. L. and Lemonte, A. J. (2018). On the Bell distribution and its
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
\code{\link{pBellBX}
}
}

\examples{
x<-rBellBX(20,2,1)
dBellBX(x,2,1)
pBellBX(x,2,1)
qBellBX(0.7,2,1)
sBellBX(x,2,1)
hBellBX(x,2,1)
mBellBX(x, 0.2,0.1, method="B")
}

