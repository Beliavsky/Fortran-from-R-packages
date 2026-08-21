\name{Random values generation from the bivariate Poisson distribution}
\alias{rbp}
\title{
Random values generation from the bivariate Poisson distribution
}

\description{
Random values generation from the bivariate Poisson distribution.
}

\usage{
rbp(n, lambda)
}

\arguments{
\item{n}{
The sample size.
}
\item{lambda}{
A vector with the three parameters, \eqn{(\lambda_1,\lambda_2,\lambda_3)} of the Poisson distribution.
}
}

\details{
In order to generate values from this distribution one needs three independent Poisson variables, \eqn{X_1 \sim \text{Po}\left(\lambda_1\right)}, \eqn{X_2 \sim \text{Po}\left(\lambda_2\right)} and \eqn{X_3 \sim \text{Po}\left(\lambda_3\right)}. Then, \eqn{\left(X,Y\right)=\left(X_1+X3,X_2+X_3\right) \sim \text{BP}\left(\lambda_1,\lambda_2,\lambda_3\right)}.
}

\value{
A matrix with n rows and 2 columns.
}

\references{
Kocherlakota S. and Kocherlakota K. (1992). Bivariate discrete distributions. CRC PRess.
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
bp.mle(x)
}

