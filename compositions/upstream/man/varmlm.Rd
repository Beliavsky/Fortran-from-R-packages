\name{varmlm}
\alias{var.lm}
\alias{var.mlm}
\title{Residual variance of a model}
\description{
  Computes the unbiased estimate for the variance of the residuals of a model.
}
\usage{
\method{var}{mlm}(x,...) 
\method{var}{lm}(x,...) 
          }
\arguments{
  \item{x}{a linear model object}
  \item{\dots}{Unused, for generic purposes only.}
}
\value{
  \item{var.lm}{returns a scalar giving the estimated variance of the residuals}
  \item{var.mlm}{returns a the estimated variance covariance matrix of  the residuals}
}

\details{
  The difference of this command to \code{var(resid(X))} is that this
  command correctly adjusts for the degrees of freedom of the model.
  
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}, Raimon
Tolosana-Delgado}

\seealso{
  \code{\link{vcov}}
}

\examples{
data(Orange)
var(lm(circumference~age,data=Orange))
var(lm(cbind(circumference,age)~age,data=Orange))

}
\keyword{multivariate}
