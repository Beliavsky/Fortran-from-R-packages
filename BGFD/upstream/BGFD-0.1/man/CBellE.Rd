\name{The complementary Bell exponential distribution}
\alias{dCBellE}
\alias{pCBellE}
\alias{qCBellE}
\alias{rCBellE}
\alias{sCBellE}
\alias{hCBellE}
\alias{mCBellE}
\title{The complementary Bell exponential distribution}

\description{Density, probability, quantile function, random generation, survival function, hazard rate function and maximum likelihood estimates from the complementary Bell exponential distribution.}

\usage{
dCBellE(x, alpha, lambda, log = FALSE)
pCBellE(x, alpha, lambda, log.p = FALSE, lower.tail = TRUE)
qCBellE(p, alpha, lambda, log.p = FALSE, lower.tail = TRUE)
rCBellE(n, alpha, lambda)
sCBellE(x, alpha, lambda, log.p = FALSE, lower.tail = TRUE)
hCBellE(x, alpha, lambda, log = FALSE,log.p = FALSE, lower.tail = TRUE)
mCBellE(x, alpha,lambda, method="B")
}

\arguments{
\item{x}{
A vector of (non-negative integer) quantiles.
}
\item{p}{
A vector of probablities.
}
\item{n}{
The number of random values to be generated under the complementary Bell exponential distribution.
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
\item{method}{the procedure for optimizing the log-likelihood function after setting the intial values of the parameters and data values for which the complementary Bell exponential distribution is fitted. It could be "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", or "SANN". "BFGS" is set as the default.
}
}

\details{The functions allow fitting the complementary Bell exponential distribution and evaluating the probability density function, cumulative distribution function, 
quantile function, random numbers, survival function, hazard rate function, and maximum likelihood estimates (MLEs) of the unknown parameters with standard error (SE) of 
the estimates. It also provides the goodness-of-fit measures such as the Akaike Information Criterion (AIC), Bayesian Information Criterion (BIC), 
the minimum value of the negative log-likelihood function, Anderson-Darling (A) test, Cramer-Von-Mises (W) test, Kolmogorov-Smirnov test, P-value and convergence status. 
}


\value{
dCBellE gives the (log) probability function.
pCBellE gives the (log) distribution function.
qCBellE gives the quantile function.
rCBellE generates random values.
sCBellE gives the survival function.
hCBellE gives the hazard rate function.
mCBellE gives the estimated parameters along with SE and goodness-of-fit measures.
}
\references{
Algarni, A. (2022). Group Acceptance Sampling Plan Based on New Compounded Three-Parameter Weibull Model. 
Axioms, 11(9): 438.

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
\code{\link{pCBellEE}
}
}

\examples{
x<-rCBellE(20,2,1)
dCBellE(x,2,1)
pCBellE(x,2,1)
qCBellE(0.7,2,1)
sCBellE(x,2,1)
hCBellE(x,2,1)
mCBellE(x, 0.2,0.1, method="B")
}

