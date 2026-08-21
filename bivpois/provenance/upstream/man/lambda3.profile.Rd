\name{Profile likelihood of the covariance parameter of the bivariate Poisson distribution}
\alias{lambda3.profile}
\title{
Profile likelihood of the covariance parameter (\eqn{\lambda_3}) of the bivariate Poisson distribution
}

\description{
Profile likelihood of the covariance parameter (\eqn{\lambda_3}) of the bivariate Poisson distribution.
}

\usage{
lambda3.profile(x1, x2 = NULL)
}

\arguments{
\item{x1}{
Either a numerical vector with the values of the first variable or a matrix with 2 columns containing both variables. In the latter case, x2 must be NULL.
}
\item{x2}{
A numerical vector with the values of the second. If x1 is a matrix with 2 columns containing both variables, x2 must be NULL.}
}

\details{
The function plots the profile log-likelihood of \eqn{\lambda_3} and computes the relevant 95\% confidence interval for the parameter.
}

\value{
A plot with the profile log-likelihood of \eqn{\lambda_3} and a vector with the 95\% confidence interval for the parameter.
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
\code{\link{bp.mle} }
}

\examples{
x <- rbp( 300, c(3, 5, 2) )
lambda3.profile(x)
}

