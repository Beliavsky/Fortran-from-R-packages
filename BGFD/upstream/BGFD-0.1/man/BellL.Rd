\name{The Bell Lomax distribution}
\alias{dBellL}
\alias{pBellL}
\alias{qBellL}
\alias{rBellL}
\alias{sBellL}
\alias{hBellL}
\alias{mBellL}
\title{The Bell Lomax distribution}

\description{Density, probability, quantile function, random generation, survival function, hazard rate function and maximum likelihood estimates from the Bell Lomax distribution.
}

\usage{
dBellL(x, b, q, lambda, log = FALSE)
pBellL(x, b, q, lambda, log.p = FALSE, lower.tail = TRUE)
qBellL(p, b, q, lambda, log.p = FALSE, lower.tail = TRUE)
rBellL(n, b, q, lambda)
sBellL(x, b, q, lambda, log.p = FALSE, lower.tail = TRUE)
hBellL(x, b, q, lambda, log = FALSE,log.p = FALSE, lower.tail = TRUE)
mBellL(x, b, q, lambda, method="B")
}

\arguments{
\item{x}{
A vector of (non-negative integer) quantiles.
}
\item{p}{
A vector of probablities.
}
\item{n}{
The number of random values to be generated under the Bell Lomax distribution.
}
\item{lambda}{
The strictly positive parameter of the Bell distribution (\eqn{\lambda > 0}).
}
\item{b}{
The strictly positive parameter of the baseline Lomax distribution (\eqn{b > 0}).
}
\item{q}{
The strictly positive shapes parameter of the baseline Lomax distribution (\eqn{q > 0}).
}
\item{lower.tail}{if FALSE then 1-F(x) are returned and quantiles are computed 1-p.
}

\item{log}{if TRUE, probabilities p are given as log(p).
}
\item{log.p}{if TRUE, probabilities p are given for exp(p).
}
\item{method}{the procedure for optimizing the log-likelihood function after setting the intial values of the parameters and data values for which the Bell Lomax distribution is fitted. It could be "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", or "SANN". "BFGS" is set as the default.
}
}

\details{The functions allow fitting the compounded Bell Lomax distribution and evaluating the probability density function, cumulative distribution function, 
quantile function, random numbers, survival function, hazard rate function, and maximum likelihood estimates (MLEs) of the unknown parameters with standard error (SE) of 
the estimates. It also provides the goodness-of-fit measures such as the Akaike Information Criterion (AIC), Bayesian Information Criterion (BIC), 
the minimum value of the negative log-likelihood function, Anderson-Darling (A) test, Cramer-Von-Mises (W) test, Kolmogorov-Smirnov test, P-value and convergence status. 
}

\value{
dBellL gives the (log) probability function.
pBellL gives the (log) distribution function.
qBellL gives the quantile function.
rBellL generates random values.
sBellL gives the survival function.
hBellL gives the hazard rate function.
mBellL gives the maximum likelihood estimates and goodness-of-fit measures.
}

\references{
Fayomi, A., Tahir, M. H., Algarni, A., Imran, M. and Jamal, F. (2022). A new useful exponential model
with applications to quality control and actuarial data. Computational Intelligence and Neuroscience.

Kleiber, C. and Kotz, S. (2003). Statistical size distributions in economics and actuarial sciences. John Wiley & Sons.

Castellares, F., Ferrari, S. L., and Lemonte, A. J. (2018). On the Bell distribution and its
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
\code{\link{pBellF}
}
}

\examples{
x<-rBellL(20,2,1,0.7)
dBellL(x,2,1,1)
pBellL(x,2,1,1)
qBellL(0.7,2,1,1)
sBellL(x,2,1,1.2)
hBellL(x,2,1,1.2)
mBellL(x, 0.2,0.1,1.2, method="B")
}

