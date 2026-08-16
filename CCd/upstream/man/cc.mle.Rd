\name{cc.mle}
\alias{cc.mle}
\alias{cc.mle0}

\title{
Maximum likelihood estimation of the CC distribution
}

\description{
Maximum likelihood estimation of the CC distribution.
}

\usage{
cc.mle(y)
cc.mle0(y, tol = 1e-7)
}

\arguments{
\item{y}{
A vector with integer values.
}
\item{tol}{
The tolerance value to terminate the maximization algorithm.
}
}

\details{
The function cc.mle0() uses the \code{\link{optimize}} function to perform MLE when the location parameter is zero, just as proposed by Papadatos (2022). The function cc.mle() uses the \code{\link{optim}} function when the location is not assumed zero.
}

\value{
A list including:
\item{param}{
For the cc.mle() a vector of the \eqn{\lambda} and \eqn{\mu} parameters.
}
\item{lambda}{
For the cc.mle0() the \eqn{\lambda} parameter.
}
\item{loglik}{
The value of the maximized log-likelihood.
}
}

\references{
Papadatos N. (2022). The characteristic function of the discrete Cauchy distribution In Memory of T. Cacoullos. Journal of Statistical Theory and Practice, 16(3): 47.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

\seealso{
\code{ \link{loc0.test}, \link{dcc}, \link{cc.reg}
}
}

\examples{
y <- round( rcauchy(100, 3, 10) )
cc.mle(y)

y <- round( rcauchy(100, 0, 10) )
cc.mle0(y)
}


