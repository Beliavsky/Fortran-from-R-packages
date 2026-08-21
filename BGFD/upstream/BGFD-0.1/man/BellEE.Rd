\name{The Bell extended exponential distribution}
\alias{dBellEE}
\alias{pBellEE}
\alias{qBellEE}
\alias{rBellEE}
\alias{sBellEE}
\alias{hBellEE}
\alias{mBellEE}
\title{The Bell extended exponential distribution
}

\description{Density, probability, quantile function, random generation, survival function, hazard rate function and maximum likelihood estimates from the Bell extended exponential distribution.}

\usage{
dBellEE(x, alpha, beta, lambda, log = FALSE)
pBellEE(x, alpha, beta, lambda, log.p = FALSE, lower.tail = TRUE)
qBellEE(p, alpha, beta, lambda, log.p = FALSE, lower.tail = TRUE)
rBellEE(n, alpha, beta, lambda)
sBellEE(x, alpha, beta, lambda, log.p = FALSE, lower.tail = TRUE)
hBellEE(x, alpha, beta, lambda, log = FALSE,log.p = FALSE, lower.tail = TRUE)
mBellEE(x, alpha, beta, lambda, method="B")
}

\arguments{
\item{x}{
A vector of (non-negative integer) quantiles.
}
\item{p}{
A vector of probablities.
}
\item{n}{
The number of random values to be generated under the Bell extended exponential distribution.
}
\item{lambda}{
The strictly positive parameter of the Bell distribution (\eqn{\lambda > 0}).
}
\item{alpha}{
The strictly positive scale parameter of the baseline extended exponential distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive shape parameter of the baseline extended exponential distribution (\eqn{\beta > 0}).
}
\item{lower.tail}{if FALSE then 1-F(x) are returned and quantiles are computed 1-p.
}

\item{log}{if TRUE, probabilities p are given as log(p).
}
\item{log.p}{if TRUE, probabilities p are given for exp(p).
}
\item{method}{the procedure for optimizing the log-likelihood function after setting the intial values of the parameters and data values for which the Bell extended exponential distribution is fitted. It could be "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", or "SANN". "BFGS" is set as the default.
}
}

\details{The functions allow fitting the compounded Bell extended exponential distribution and evaluating the probability density function, cumulative distribution function, 
quantile function, random numbers, survival function, hazard rate function, and maximum likelihood estimates (MLEs) of the unknown parameters with standard error (SE) of 
the estimates. It also provides the goodness-of-fit measures such as the Akaike Information Criterion (AIC), Bayesian Information Criterion (BIC), 
the minimum value of the negative log-likelihood function, Anderson-Darling (A) test, Cramer-Von-Mises (W) test, Kolmogorov-Smirnov test, P-value and convergence status. 
}


\value{
dBellEE gives the (log) probability function.
pBellEE gives the (log) distribution function.
qBellEE gives the quantile function.
rBellEE generates random values.
sBellEE gives the survival function.
hBellEE gives the hazard rate function.
mBellEE gives the estimated parameters along with SE and goodness-of-fit measures.
}

\references{
Fayomi, A., Tahir, M. H., Algarni, A., Imran, M. and Jamal, F. (2022). A new useful exponential model
with applications to quality control and actuarial data.
Computational Intelligence and Neuroscience, 2022.

Nadarajah, S. (2011). The exponentiated exponential distribution: a survey. 
AStA Advances in Statistical Analysis, 95: 219-251.

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
\code{\link{pBellE}
}
}

\examples{
x<-rBellEE(20,2,1,1)
dBellEE(x,2,1,1.2)
pBellEE(x,2,1,1.2)
qBellEE(0.7,2,1,1.2)
sBellEE(x,2,1,1.2)
hBellEE(x,2,1,1.2)
mBellEE(x, 0.2,0.1,1.2, method="B")
}

