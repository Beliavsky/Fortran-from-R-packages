\name{The complementary Bell Burr-X distribution}
\alias{dCBellBX}
\alias{pCBellBX}
\alias{qCBellBX}
\alias{rCBellBX}
\alias{sCBellBX}
\alias{hCBellBX}
\alias{mCBellBX}
\title{The complementary Bell Burr-X distribution}

\description{Density, probability, quantile function, random generation, survival function, hazard rate function and maximum likelihood estimates from the complementary Bell Burr-X distribution.
}

\usage{
dCBellBX(x, a, lambda, log = FALSE)
pCBellBX(x, a, lambda, log.p = FALSE, lower.tail = TRUE)
qCBellBX(p, a, lambda, log.p = FALSE, lower.tail = TRUE)
rCBellBX(n, a, lambda)
sCBellBX(x, a, lambda, log.p = FALSE, lower.tail = TRUE)
hCBellBX(x, a, lambda, log = FALSE, log.p = FALSE, lower.tail = TRUE)
mCBellBX(x, a, lambda, method="B")
}

\arguments{
\item{x}{
A vector of (non-negative integer) quantiles.
}
\item{p}{
A vector of probablities.
}
\item{n}{
The number of random values to be generated under the complementary Bell Burr-X distribution.
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
\item{method}{the procedure for optimizing the log-likelihood function after setting the intial values of the parameters and data values for which the complementary Bell Burr-X distribution is fitted. It could be "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", or "SANN". "BFGS" is set as the default.
}
}

\details{The functions allow fitting the complementary Bell Burr-X distribution and evaluating the probability density function, cumulative distribution function, 
quantile function, random numbers, survival function, hazard rate function, and maximum likelihood estimates (MLEs) of the unknown parameters with standard error (SE) of 
the estimates. It also provides the goodness-of-fit measures such as the Akaike Information Criterion (AIC), Bayesian Information Criterion (BIC), 
the minimum value of the negative log-likelihood function, Anderson-Darling (A) test, Cramer-Von-Mises (W) test, Kolmogorov-Smirnov test, P-value and convergence status. 
}

\value{
dCBellBX gives the (log) probability function.
pCBellBX gives the (log) distribution function.
qCBellBX gives the quantile function.
rCBellBX generates random values.
sCBellBX gives the survival function.
hCBellBX gives the hazard rate function.
mCBellBX gives the estimated parameters along with SE and goodness-of-fit measures.
}

\references{
Algarni, A. (2022). Group Acceptance Sampling Plan Based on New Compounded Three-Parameter Weibull Model. 
Axioms, 11(9): 438.

Zimmer, W. J., Keats, J. B. and Wang, F. K. (1998). The Burr XII distribution in reliability analysis. 
Journal of Quality Technology, 30(4), 386-394.

Castellares, F., Ferrari, S. L. and Lemonte, A. J. (2018). On the bell distribution and its
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
\code{\link{pBellB}
}
}

\examples{
x<-rCBellBX(20,0.2,1)
dCBellBX(x,2,1)
pCBellBX(x,2,1)
qCBellBX(0.7,2,1)
sCBellBX(x,2,1)
hCBellBX(x,2,1)
mCBellBX(x, 0.2,0.1, method="B")
}

