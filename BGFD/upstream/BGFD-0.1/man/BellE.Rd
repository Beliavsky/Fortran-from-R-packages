\name{The Bell exponential distribution}
\alias{dBellE}
\alias{pBellE}
\alias{qBellE}
\alias{rBellE}
\alias{sBellE}
\alias{hBellE}
\alias{mBellE}
\title{The Bell exponential distribution}

\description{Density, probability, quantile function, random generation, survival function, hazard rate function and maximum likelihood estimates from the Bell exponential distribution.}

\usage{
dBellE(x, alpha, lambda, log = FALSE)
pBellE(x, alpha, lambda, log.p = FALSE, lower.tail = TRUE)
qBellE(p, alpha, lambda, log.p = FALSE, lower.tail = TRUE)
rBellE(n, alpha, lambda)
sBellE(x, alpha, lambda, log.p = FALSE, lower.tail = TRUE)
hBellE(x, alpha, lambda, log = FALSE,log.p = FALSE, lower.tail = TRUE)
mBellE(x, alpha,lambda, method="B")
}

\arguments{
\item{x}{
A vector of (non-negative integer) quantiles.
}
\item{p}{
A vector of probablities.
}
\item{n}{
The number of random values to be generated under the Bell exponential distribution.
}
\item{lambda}{
The strictly positive parameter of the Bell distribution (\eqn{\lambda > 0}).
}
\item{alpha}{
The strictly positive parameter of the baseline exponential distribution (\eqn{\alpha > 0}).
}
\item{lower.tail}{if FALSE then 1-F(x) are returned and quantiles are computed 1-p.
}

\item{log}{if TRUE, probabilities p are given as log(p).
}
\item{log.p}{if TRUE, probabilities p are given for exp(p).
}
\item{method}{the procedure for optimizing the log-likelihood function after setting the intial values of the parameters and data values for which the Bell exponential distribution is fitted. It could be "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", or "SANN". "BFGS" is set as the default.
}
}

\details{The functions allow fitting the compounded Bell exponential distribution and evaluating the probability density function, cumulative distribution function, 
quantile function, random numbers, survival function, hazard rate function, and maximum likelihood estimates (MLEs) of the unknown parameters with standard error (SE) of 
the estimates. It also provides the goodness-of-fit measures such as the Akaike Information Criterion (AIC), Bayesian Information Criterion (BIC), 
the minimum value of the negative log-likelihood function, Anderson-Darling (A) test, Cramer-Von-Mises (W) test, Kolmogorov-Smirnov test, P-value and convergence status. 
}

\value{
dBellE gives the (log) probability function.
pBellE gives the (log) distribution function.
qBellE gives the quantile function.
rBellE generates random values.
sBellE gives the survival function.
hBellE gives the hazard rate function.
mBellE gives the estimated parameters along with SE and goodness-of-fit measures.
}

\references{
Fayomi, A., Tahir, M. H., Algarni, A., Imran, M. and Jamal, F. (2022). A new useful exponential model
with applications to quality control and actuarial data.
Computational Intelligence and Neuroscience.

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
\code{\link{pBellW}
}
}

\examples{
x<-rBellE(20,2,1)
dBellE(x,2,1)
pBellE(x,2,1)
qBellE(0.7,2,1)
sBellE(x,2,1)
hBellE(x,2,1)
mBellE(x, 0.2,0.1, method="B")
}

