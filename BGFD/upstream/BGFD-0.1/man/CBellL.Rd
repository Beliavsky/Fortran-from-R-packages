\name{The complementary Bell Lomax  distribution}
\alias{dCBellL}
\alias{pCBellL}
\alias{qCBellL}
\alias{rCBellL}
\alias{sCBellL}
\alias{hCBellL}
\alias{mCBellL}
\title{The complementary Bell Lomax  distribution}

\description{Density, probability, quantile function, random generation, survival function, hazard rate function and maximum likelihood estimates from the complementary Bell Lomax  distribution.}

\usage{
dCBellL(x, b, q, lambda, log = FALSE)
pCBellL(x, b, q, lambda, log.p = FALSE, lower.tail = TRUE)
qCBellL(p, b, q, lambda, log.p = FALSE, lower.tail = TRUE)
rCBellL(n, b, q, lambda)
sCBellL(x, b, q, lambda, log.p = FALSE, lower.tail = TRUE)
hCBellL(x, b, q, lambda, log = FALSE,log.p = FALSE, lower.tail = TRUE)
mCBellL(x, b, q, lambda, method="B")
}

\arguments{
\item{x}{
A vector of (non-negative integer) quantiles.
}
\item{p}{
A vector of probablities.
}
\item{n}{
The number of random values to be generated under the complementary Bell Lomax  distribution.
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

\item{log}{if TRUE, probabilities p are given as log(p)
}
\item{log.p}{if TRUE, probabilities p are given for exp(p)
}
\item{method}{the procedure for optimizing the log-likelihood function after setting the intial values of the parameters and data values for which the complementary Bell Lomax distribution is fitted. It could be "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", or "SANN". "BFGS" is set as the default.
}
}

\details{The functions allow fitting the complementary Bell Lomax distribution and evaluating the probability density function, cumulative distribution function, 
quantile function, random numbers, survival function, hazard rate function, and maximum likelihood estimates (MLEs) of the unknown parameters with standard error (SE) of 
the estimates. It also provides the goodness-of-fit measures such as the Akaike Information Criterion (AIC), Bayesian Information Criterion (BIC), 
the minimum value of the negative log-likelihood function, Anderson-Darling (A) test, Cramer-Von-Mises (W) test, Kolmogorov-Smirnov test, P-value and convergence status. 
}

\value{
dCBellL gives the (log) probability function.
pCBellL gives the (log) distribution function.
qCBellL gives the quantile function.
rCBellL generates random values.
sCBellL gives the survival function.
hCBellL gives the hazard rate function.
mCBellL gives the maximum likelihood estimates and goodness-of-fit measures.
}

\references{
Algarni, A. (2022). Group Acceptance Sampling Plan Based on New Compounded Three-Parameter Weibull Model. 
Axioms, 11(9): 438.

Kleiber, C. and Kotz, S. (2003). Statistical size distributions in economics and actuarial sciences. 
John Wiley & Sons.

Castellares, F., Ferrari, S. L. and Lemonte, A. J. (2018). On the Bell distribution and its
associated regression model for count data. Applied Mathematical Modelling, 56, 172-185.
}


\author{
Muhammad Imran and Michail Tsagris.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com} and Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{pCBellF}
}
}

\examples{
x<-rCBellL(20,2,1,1.2)
dCBellL(x,2,1,0.5)
pCBellL(x,2,1,0.5)
qCBellL(0.7,2,1,0.5)
sCBellL(x,2,1,0.5)
hCBellL(x,2,1,0.5)
mCBellL(x, 0.2,0.1,0.8, method="B")
}

