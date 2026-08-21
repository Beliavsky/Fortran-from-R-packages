\name{skewness}
\alias{skewness}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Skewness}
\description{
Computes the skewness of the data} 
\usage{
skewness(x, na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{data}
  \item{na.rm}{logical variable, if true, missing values are excluded from analysis}
}
\value{The function returns the skewness of the data.}
\details{Given data \eqn{x_1,x_2,\ldots, x_n}, the sample skewness is 
defined by the formula:
\deqn{\frac{\sum_{i=1}^n (x_i-\bar{x})^3/n}{(\sum_{i=1}^n 
(x_i-\bar{x})^2/n)^{3/2}}.}
}
\author{Kung-Sik Chan}
\examples{
data(CREF)
r.cref=diff(log(CREF))*100
skewness(r.cref)
}
\keyword{methods}
