\name{cc.reg}
\alias{cc.reg}

\title{
Regression modelling with the CC distribution
}

\description{
Regression modelling with the CC distribution.
}

\usage{
cc.reg(y, x, tol = 1e-6)
}

\arguments{
\item{y}{
The response variable, a vector with integer values.
}
\item{x}{
A vector or matrix with with the predictor variables.
}
\item{tol}{
The tolerance value to terminate the maximization algorithm.
}
}

\details{
Regression modelling assuming that the counts follow the CC distribution is implemented.
}

\value{
A list including:
\item{lambda}{
The \eqn{\lambda} parameter.
}
\item{be}{
The regression coefficients.
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
\code{ \link{cc.mle}
}
}

\examples{
y <- round( rcauchy(150, 3, 10) )
x <- iris[, 1:2]
cc.reg(y, x)
}


