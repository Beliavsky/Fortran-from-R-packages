\name{The complementary Bell Burr-12 distribution}
\alias{dCBellB}
\alias{pCBellB}
\alias{qCBellB}
\alias{rCBellB}
\alias{sCBellB}
\alias{hCBellB}
\alias{mCBellB}
\title{The complementary Bell Burr-12 distribution}

\description{Density, probability, quantile function, random generation, survival function, hazard rate function and maximum likelihood estimates from the complementary Bell Burr-12 distribution.
}

\usage{
dCBellB(x, a, b, k, lambda, log = FALSE)
pCBellB(x, a, b, k, lambda, log.p = FALSE, lower.tail = TRUE)
qCBellB(p, a, b, k, lambda, log.p = FALSE, lower.tail = TRUE)
rCBellB(n, a, b, k, lambda)
sCBellB(x, a, b, k, lambda, log.p = FALSE, lower.tail = TRUE)
hCBellB(x, a, b, k, lambda, log = FALSE,log.p = FALSE, lower.tail = TRUE)
mCBellB(x, a, b, k, lambda, method="B")
}

\arguments{
\item{x}{
A vector of (non-negative integer) quantiles.
}
\item{p}{
A vector of probablities.
}
\item{n}{
The number of random values to be generated under the complementary Bell Burr-12 distribution.
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
\item{lower.tail}{if FALSE then 1-F(x) are returned and quantiles are computed 1-p.
}

\item{log}{if TRUE, probabilities p are given as log(p).
}
\item{log.p}{if TRUE, probabilities p are given for exp(p).
}
\item{method}{the procedure for optimizing the log-likelihood function after setting the intial values of the parameters and data values for which the complementary Bell Burr-12 distribution is fitted. It could be "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", or "SANN". "BFGS" is set as the default.
}
}

\details{The functions allow fitting the complementary Bell Burr-12 distribution and evaluating the probability density function, cumulative distribution function, 
quantile function, random numbers, survival function, hazard rate function, and maximum likelihood estimates (MLEs) of the unknown parameters with standard error (SE) of 
the estimates. It also provides the goodness-of-fit measures such as the Akaike Information Criterion (AIC), Bayesian Information Criterion (BIC), 
the minimum value of the negative log-likelihood function, Anderson-Darling (A) test, Cramer-Von-Mises (W) test, Kolmogorov-Smirnov test, P-value and convergence status. 
}

\value{
dCBellB gives the (log) probability function.
pCBellB gives the (log) distribution function.
qCBellB gives the quantile function.
rCBellB generates random values.
sCBellB gives the survival function.
hCBellB gives the hazard rate function.
mCBellB gives the estimated parameters along with SE and goodness-of-fit measures.
}

\references{
Algarni, A. (2022). Group Acceptance Sampling Plan Based on New Compounded Three-Parameter Weibull Model. 
Axioms, 11(9): 438.

Zimmer, W. J., Keats, J. B. and Wang, F. K. (1998). The Burr XII distribution in reliability analysis. 
Journal of Quality Technology, 30(4), 386-394.

Castellares, F., Ferrari, S. L. and Lemonte, A. J. (2018). On the Bell distribution and its
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
\code{\link{pBellBX}
}
}

\examples{
x<-rCBellB(20,2,1,0.7,1)
dCBellB(x,2,1,1,0.2)
pCBellB(x,2,1,1,0.2,1)
qCBellB(0.7,2,1,1,0.2,1)
sCBellB(x,2,1,1,0.2,1)
hCBellB(x,2,1,1,0.2,1)
mCBellB(x, 0.2,0.1,1.2,0.2,method="B")
}

