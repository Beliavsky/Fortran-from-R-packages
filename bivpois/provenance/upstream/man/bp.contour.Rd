\name{Contour plot of the bivariate Poisson distribution}
\alias{bp.contour}
\title{
Contour plot of the bivariate Poisson distribution
}

\description{
Contour plot of the bivariate Poisson distribution.
}

\usage{
bp.contour(x1, x2 = NULL, lambda)
}

\arguments{
\item{x1}{
Either a numerical vector with the values of the first variable or a matrix with 2 columns containing both variables. In the latter case, x2 must be NULL.
}
\item{x2}{
A numerical vector with the values of the second. If x1 is a matrix with 2 columns containing both variables, x2 must be NULL.}
\item{lambda}{
A vector with three numbers, the estimates of the \eqn{\lambda_s}.
}
}

\details{
This is actually a scatter plot with the probability values computed at each pair of points.
}

\value{
The contour plot of the bivariate Poisson distribution.
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
lambda <- bp.mle(x)$lambda
bp.contour(x, lambda = lambda)
}

