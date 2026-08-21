\name{Density computation of the bivariate Poisson distribution}
\alias{dbp}
\title{
Density computation of the bivariate Poisson distribution
}

\description{
Density computation of the bivariate Poisson distribution.
}

\usage{
dbp(x1, x2 = NULL, lambda, logged = TRUE)
}

\arguments{
\item{x1}{
Either a numerical vector with the values of the first variable or a matrix with 2 columns containing both variables. In the latter case, x2 must be NULL.
}
\item{x2}{
A numerical vector with the values of the second. If x1 is a matrix with 2 columns containing both variables, x2 must be NULL.
}
\item{lambda}{
A vector with three numbers, the estimates of the \eqn{\lambda_s}.
}
\item{logged}{
Should the logarithm of the density values be computed? The default value is TRUE.
}
}

\details{
The density of the bivariate Poisson distribution is computed.
}

\value{
A vector with the logged density values.
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
a <- bp.mle(x)
f <- dbp(x, lambda = a$lambda)
sum(f)
}

