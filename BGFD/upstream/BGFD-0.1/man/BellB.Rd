\name{The Bell Burr-12 distribution}

\alias{dBellB}
\alias{pBellB}
\alias{qBellB}
\alias{rBellB}
\alias{sBellB}
\alias{hBellB}
\alias{mBellB}

\title{The Bell Burr-12 distribution
}

\description{Density, probability, quantile function, random generation, survival function, hazard rate function and maximum likelihood estimates from the Bell Burr-12 distribution.
}

\usage{
dBellB(x, a, b, k, lambda, log = FALSE)
pBellB(x, a, b, k, lambda, log.p = FALSE, lower.tail = TRUE)
qBellB(p, a, b, k, lambda, log.p = FALSE, lower.tail = TRUE)
rBellB(n, a, b, k, lambda)
sBellB(x, a, b, k, lambda, log.p = FALSE, lower.tail = TRUE)
hBellB(x, a, b, k, lambda, log = FALSE,log.p = FALSE, lower.tail = TRUE)
mBellB(x, a, b, k, lambda, method="B")
}

\arguments{
\item{x}{
A vector of (non-negative integer) quantiles.
}
\item{p}{
A vector of probablities.
}
\item{n}{
The number of random values to be generated under the Bell Burr-12 distribution.
}
\item{lambda}{
The strictly positive parameter of the Bell distribution (\eqn{\lambda > 0}).
}
\item{a}{
The strictly positive scale parameter of the baseline Burr-12 distribution (\eqn{a > 0}).
}
\item{b}{
The strictly positive shape parameter of the baseline Burr-12 distribution (\eqn{b > 0}).
}
\item{k}{
The strictly positive shape parameter of the baseline Burr-12 distribution (\eqn{k > 0}).
}
\item{lower.tail}{
If FALSE then \eqn{1-F(x)} are returned and quantiles are computed \eqn{1-p}.
}
\item{log}{
If TRUE, probabilities p are given as log(p).
}
\item{log.p}{
If TRUE, probabilities p are given for exp(p).
}
\item{method}{the procedure for optimizing the log-likelihood function after setting the intial values of the parameters and data values for which the Bell Burr-12 distribution is fitted. It could be "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", or "SANN". "BFGS" is set as the default.
}
}

\details{The functions allow fitting the compounded Bell-Burr-12 distribution and evaluating the probability density function, cumulative distribution function,
quantile function, random numbers, survival function, hazard rate function, and maximum likelihood estimates (MLEs) of the unknown parameters with standard error (SE) of
the estimates. It also provides the goodness-of-fit measures such as the Akaike Information Criterion (AIC), Bayesian Information Criterion (BIC),
the minimum value of the negative log-likelihood function, Anderson-Darling (A) test, Cramer-Von-Mises (W) test, Kolmogorov-Smirnov test, P-value and convergence status.
}

\value{
dBellB gives the (log) probability function.
pBellB gives the (log) distribution function.
qBellB gives the quantile function.
rBellB generates random values.
sBellB gives the survival function.
hBellB gives the hazard rate function.
mBellB gives the estimated parameters along with SE and goodness-of-fit measures.
}

\references{
Fayomi, A., Tahir, M. H., Algarni, A., Imran, M. and Jamal, F. (2022). A new useful exponential model
with applications to quality control and actuarial data.
Computational Intelligence and Neuroscience.

Zimmer, W. J., Keats, J. B. and Wang, F. K. (1998). The Burr XII distribution in reliability analysis.
Journal of Quality Technology, 30(4): 386-394.

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
\code{\link{pBellEW}
}
}

\examples{
x<-rBellB(20,2,0.4,1.2,1.5)
dBellB(x,2,1,2,1.5)
pBellB(x,2,1,2,1.5)
qBellB(0.7,2,1,2,1.5)
sBellB(x,2,1,2,1.5)
hBellB(x,2,1,2,1.5)
mBellB(x, 0.2,0.4,1.5,1.2, method="B")
}

