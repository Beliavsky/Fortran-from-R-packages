\name{kurtosis}
\alias{kurtosis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Kurtosis}
\description{
Computes the Kurtosis.
}
\usage{
kurtosis(x, na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{data}
  \item{na.rm}{logical variable, if true, missing values are excluded from analysis}
}
\value{The function returns the kurtosis of the data.}
\details{Given data \eqn{x_1,x_2,\ldots, x_n}, the sample kurtosis is 
defined by the formula:
\deqn{\frac{\sum_{i=1}^n (x_i-\bar{x})^4/n}{(\sum_{i=1}^n 
(x_i-\bar{x})^2/n)^2}-3.}
}
\author{Kung-Sik Chan}
\examples{
data(CREF)
r.cref=diff(log(CREF))*100
kurtosis(r.cref)
}
\keyword{methods}
