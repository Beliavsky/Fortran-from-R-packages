\name{The complementary Bell exponentiated Weibull distribution}
\alias{dCBellEW}
\alias{pCBellEW}
\alias{qCBellEW}
\alias{rCBellEW}
\alias{sCBellEW}
\alias{hCBellEW}
\alias{mCBellEW}
\title{The complementary Bell exponentiated Weibull distribution}

\description{Density, probability, quantile function, random generation, survival function, hazard rate function and maximum likelihood estimates from the complementary Bell exponentiated Weibull distribution.}
\usage{
dCBellEW(x, alpha, beta, theta, lambda, log = FALSE)
pCBellEW(x, alpha, beta, theta, lambda, log.p = FALSE, lower.tail = TRUE)
qCBellEW(p, alpha, beta, theta, lambda, log.p = FALSE, lower.tail = TRUE)
rCBellEW(n, alpha, beta, theta, lambda)
sCBellEW(x, alpha, beta, theta, lambda, log.p = FALSE, lower.tail = TRUE)
hCBellEW(x, alpha, beta, theta, lambda, log = FALSE,log.p = FALSE, lower.tail = TRUE)
mCBellEW(x, alpha, beta, theta, lambda, method="B")
}

\arguments{
\item{x}{
A vector of (non-negative integer) quantiles.
}
\item{p}{
A vector of probablities.
}
\item{n}{
The number of random values to be generated under the complementary Bell exponentiated Weibull distribution.
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
\item{method}{the procedure for optimizing the log-likelihood function after setting the intial values of the parameters and data values for which the complementary Bell exponentiated Weibull distribution is fitted. It could be "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", or "SANN". "BFGS" is set as the default.
}
}

\details{The functions allow fitting the complementary Bell exponentiated Weibull distribution and evaluating the probability density function, cumulative distribution function, 
quantile function, random numbers, survival function, hazard rate function, and maximum likelihood estimates (MLEs) of the unknown parameters with standard error (SE) of 
the estimates. It also provides the goodness-of-fit measures such as the Akaike Information Criterion (AIC), Bayesian Information Criterion (BIC), 
the minimum value of the negative log-likelihood function, Anderson-Darling (A) test, Cramer-Von-Mises (W) test, Kolmogorov-Smirnov test, P-value and convergence status. 
}

\value{
dCBellEW gives the (log) probability function.
pCBellEW gives the (log) distribution function.
qCBellEW gives the quantile function.
rCBellEW generates random values.
sCBellEW gives the survival function.
hCBellEW gives the hazard rate function.
mCBellEW gives the estimated parameters along with SE and goodness-of-fit measures.
}

\references{
Algarni, A. (2022). Group Acceptance Sampling Plan Based on New Compounded Three-Parameter Weibull Model. 
Axioms, 11(9): 438.

Nadarajah, S., Cordeiro, G. M. and Ortega, E. M. (2013). The exponentiated Weibull distribution: a survey. 
Statistical Papers, 54, 839-877.

Castellares, F., Ferrari, S. L., and Lemonte, A. J. (2018). On the bell distribution and its
associated regression model for count data. 
Applied Mathematical Modelling, 56, 172-185.
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
x<-rCBellEW(20,2,1,0.2,2.2)
dCBellEW(x,2,1,0.5,0.2)
pCBellEW(x,2,1,0.5,0.2)
qCBellEW(0.7,2,1,0.5,0.1)
sCBellEW(x,2,1,0.5,0.2)
hCBellEW(x,2,1,0.5,0.5)
mCBellEW(x, 0.2,0.1,0.8,0.5, method="B")
}

