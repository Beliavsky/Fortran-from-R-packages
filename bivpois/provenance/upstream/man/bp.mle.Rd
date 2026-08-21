\name{Maximum likelihood estimation of the bivariate Poisson distribution}
\alias{bp.mle}
\alias{bp.mle2}
\title{
Maximum likelihood estimation of the bivariate Poisson distribution
}

\description{
Maximum likelihood estimation of the bivariate Poisson distribution.
}

\usage{
bp.mle(x1, x2 = NULL)
bp.mle2(x1, x2 = NULL)
}

\arguments{
\item{x1}{
Either a numerical vector with the values of the first variable or a matrix with 2 columns containing both variables. In the latter case, x2 must be NULL.
}
\item{x2}{
A numerical vector with the values of the second. If x1 is a matrix with 2 columns containing both variables, x2 must be NULL.}
}

\details{
Using the addition method (see function \code{\link{rbp}}) to simulate random values from the bivariate Poisson, its representation is given by
\deqn{
P\left(X=x,Y=y \right)e^{-\left(\lambda_1+\lambda_2+\lambda_3\right)}\frac{\lambda_1^x}{x!}\frac{\lambda_2^y}{y!}
\sum_{k=0}^{\min{\left(x,y\right)}} \binom {x} {k} \binom {y} {k} k! \left(\frac{\lambda_3}{\lambda_1\lambda_2}\right)^k.
}

The above form is found in Karlis and Ntzoufras (2003). This bivariate distribution allows for dependence between the two random variables. Marginally each random variable follows a Poisson distribution with \eqn{E\left(X\right)=\lambda_1+\lambda_3} and \eqn{E\left(Y\right)=\lambda_2+\lambda_3}. In addition, \eqn{Cov\left(X,Y\right)=\lambda_3}. If \eqn{\lambda_3=0}, the above expression becomes a product of two Poisson distributions. Hence, \eqn{\lambda_3} is a measure of dependence between the two random variables.

The function bp.mle() returns a lot of information and is slower than bp.mle2(), which returns fewer information, but is faster.
}

\value{
For the function bp.mle() a list including:
\item{lambda}{
A vector with the estimated values of \eqn{(\lambda_1}, \eqn{\lambda_2)} and \eqn{\lambda_3}.
Note that \eqn{\hat{\lambda}_1=\bar{x}_1 - \lambda_3} and \eqn{\hat{\lambda}_1=\bar{x}_1 - \lambda_3}, where
\eqn{\bar{x}_1} and \eqn{\bar{x}_2} are the two sample means.
}
\item{rho}{
The estimated correlation coefficient, that is: \eqn{\dfrac{\hat{\lambda}_3}{\sqrt{\left(\hat{\lambda}_1 + \hat{\lambda_3}\right)\left(\hat{\lambda}_2 + \hat{\lambda_3}\right)}}}.
}
\item{ci}{
The 95\% Confidence intervals using the observed and the asymptotic information matrix.
}
\item{loglik}{
The log-likelihood values assuming independence (\eqn{\lambda_3=0}) and assuming the bivariate Poisson distribution.
}
\item{pvalue}{
Three p-values for testing \eqn{\lambda_3=0}. These are based on the log-likelihood ratio and two Wald tests using the observed and the asymptotic information matrix.
}

For the function bp.mle2() a list including:
\item{lambda}{
A vector with the estimated values of \eqn{(\lambda_1}, \eqn{\lambda_2)} and \eqn{\lambda_3}.
Note that \eqn{\hat{\lambda}_1=\bar{x}_1 - \lambda_3} and \eqn{\hat{\lambda}_1=\bar{x}_1 - \lambda_3}, where
\eqn{\bar{x}_1} and \eqn{\bar{x}_2} are the two sample means.
}
\item{loglik}{
The log-likelihood values assuming independence (\eqn{\lambda_3=0}) and assuming the bivariate Poisson distribution.
}
}

\references{
Kawamura K. (1984). Direct calculation of maximum likelihood estimator for the bivariate
Poisson distribution. Kodai Mathematical Journal, 7(2): 211--221.

Kocherlakota S. and Kocherlakota K. (1992). Bivariate discrete distributions. CRC Press.

Karlis D. and Ntzoufras I. (2003). Analysis of sports data by using bivariate poisson models.
Journal of the Royal Statistical Society: Series D (The Statistician), 52(3): 381--393.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{rbp} }
}

\examples{
x <- rbp( 300, c(3, 5, 2) )
bp.mle(x)
}

