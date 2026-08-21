\name{The complementary Bell Fisk distribution}
\alias{dCBellF}
\alias{pCBellF}
\alias{qCBellF}
\alias{rCBellF}
\alias{sCBellF}
\alias{hCBellF}
\alias{mCBellF}
\title{The complementary Bell Fisk distribution
}

\description{Density, probability, quantile function, random generation, survival function, hazard rate function and maximum likelihood estimates from the complementary Bell Fisk distribution.}
\usage{
dCBellF(x, a, b, lambda, log = FALSE)
pCBellF(x, a, b, lambda, log.p = FALSE, lower.tail = TRUE)
qCBellF(p, a, b, lambda, log.p = FALSE, lower.tail = TRUE)
rCBellF(n, a, b, lambda)
sCBellF(x, a, b, lambda, log.p = FALSE, lower.tail = TRUE)
hCBellF(x, a, b, lambda, log = FALSE,log.p = FALSE, lower.tail = TRUE)
mCBellF(x, a, b, lambda, method="B")
}
\arguments{
\item{x}{
A vector of (non-negative integer) quantiles.
}
\item{p}{
A vector of probablities.
}
\item{n}{
The number of random values to be generated under the complementary Bell Fisk distribution.
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
\item{method}{the procedure for optimizing the log-likelihood function after setting the intial values of the parameters and data values for which the complementary Bell Fisk distribution is fitted. It could be "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", or "SANN". "BFGS" is set as the default.
}
}

\details{The functions allow fitting the complementary Bell Fisk distribution and evaluating the probability density function, cumulative distribution function, 
quantile function, random numbers, survival function, hazard rate function, and maximum likelihood estimates (MLEs) of the unknown parameters with standard error (SE) of 
the estimates. It also provides the goodness-of-fit measures such as the Akaike Information Criterion (AIC), Bayesian Information Criterion (BIC), 
the minimum value of the negative log-likelihood function, Anderson-Darling (A) test, Cramer-Von-Mises (W) test, Kolmogorov-Smirnov test, P-value and convergence status. 
}

\value{
dCBellF gives the (log) probability function.
pCBellF gives the (log) distribution function.
qCBellF gives the quantile function.
rCBellF generates random values.
sCBellF gives the survival function.
hCBellF gives the hazard rate function.
mCBellF gives the estimated parameters along with SE and goodness-of-fit measures.
}

\references{
Algarni, A. (2022). Group Acceptance Sampling Plan Based on New Compounded Three-Parameter Weibull Model. 
Axioms, 11(9): 438.

Kleiber, C. and Kotz, S. (2003). Statistical size distributions in economics and actuarial sciences. 
John Wiley & Sons.

Castellares, F., Ferrari, S. L., and Lemonte, A. J. (2018). On the Bell distribution and its
associated regression model for count data. Applied Mathematical Modelling, 56, 172- 185.
}

\author{
Muhammad Imran and Michail Tsagris.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com} and Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{pBellL}
}
}

\examples{
x<-rCBellF(20,2,1,0.2)
dCBellF(x,2,1,0.5)
pCBellF(x,2,1,0.5)
qCBellF(0.7,2,1,0.5)
sCBellF(x,2,1,0.5)
hCBellF(x,2,1,0.5)
mCBellF(x, 0.2,0.1,0.8, method="B")
}

