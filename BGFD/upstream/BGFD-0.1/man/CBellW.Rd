\name{The complementary Bell Weibull distribution}
\alias{dCBellW}
\alias{pCBellW}
\alias{qCBellW}
\alias{rCBellW}
\alias{sCBellW}
\alias{hCBellW}
\alias{mCBellW}
\title{The complementary Bell Weibull distribution}

\description{Density, probability, quantile function, random generation, survival function, hazard rate function and maximum likelihood estimates from the complementary Bell Weibull distribution.
}

\usage{
dCBellW(x, alpha, beta, lambda, log = FALSE)
pCBellW(x, alpha, beta, lambda, log.p = FALSE, lower.tail = TRUE)
qCBellW(p, alpha, beta, lambda, log.p = FALSE, lower.tail = TRUE)
rCBellW(n, alpha, beta, lambda)
sCBellW(x, alpha, beta, lambda, log.p = FALSE, lower.tail = TRUE)
hCBellW(x, alpha, beta, lambda, log = FALSE,log.p = FALSE, lower.tail = TRUE)
mCBellW(x, alpha, beta, lambda, method="B")
}

\arguments{
\item{x}{
A vector of (non-negative integer) quantiles.
}
\item{p}{
A vector of probablities.
}
\item{n}{
The number of random values to be generated under the complementary Bell Weibull distribution.
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
\item{method}{the procedure for optimizing the log-likelihood function after setting the intial values of the parameters and data values for which the complementary Bell Weibull  distribution is fitted. It could be "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", or "SANN". "BFGS" is set as the default.
}
}

\details{The functions allow fitting the complementary Bell Weibull distribution and evaluating the probability density function, cumulative distribution function, 
quantile function, random numbers, survival function, hazard rate function, and maximum likelihood estimates (MLEs) of the unknown parameters with standard error (SE) of 
the estimates. It also provides the goodness-of-fit measures such as the Akaike Information Criterion (AIC), Bayesian Information Criterion (BIC), 
the minimum value of the negative log-likelihood function, Anderson-Darling (A) test, Cramer-Von-Mises (W) test, Kolmogorov-Smirnov test, P-value and convergence status. 
}

\value{
dCBellW gives the (log) probability function.
pCBellW gives the (log) distribution function.
qCBellW gives the quantile function.
rCBellW generates random values.
sCBellW gives the survival function.
hCBellW gives the hazard rate function.
mCBellW gives the estimated parameters along with SE and goodness-of-fit measures.
}

\references{

Algarni, A. (2022). Group Acceptance Sampling Plan Based on New Compounded Three-Parameter Weibull Model. 
Axioms, 11(9): 438.

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
\code{\link{pCBellEW}
}
}

\examples{
x<-rCBellW(20,2,1,0.2)
dCBellW(x,2,1,0.5)
pCBellW(x,2,1,0.5)
qCBellW(0.7,2,1,0.5)
sCBellW(x,2,1,0.5)
hCBellW(x,2,1,0.5)
mCBellW(x, 0.2,0.1,0.8, method="B")
}

