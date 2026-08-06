\name{change.origin}
\alias{change.origin}
\title{Change Origin for a Polynomial}
\description{
  Calculate the coefficients of a polynomial relative to a new origin
  on the x axis.
}
\usage{
change.origin(p, o)
}
\arguments{
  \item{p}{an object of class \code{"polynomial"}.}
  \item{o}{a numeric scalar representing the new origin on the original
    scale.}
}
\details{
  Let \eqn{P(x) = \sum_i p_i x^i} be a given polynomial and consider
  writing \eqn{P(x) = \sum_j q_j (x - o)^j}.  This function calculates
  the coefficients \eqn{q_j} and returns the result as a polynomial.
}
\value{
  A polynomial with coefficients relative to the re-located x axis.
}
\examples{
pr <- poly.calc(1:5)
pr
## -120 + 274*x - 225*x^2 + 85*x^3 - 15*x^4 + x^5
change.origin(pr, 3)
## 4*x - 5*x^3 + x^5
}
\keyword{symbolmath}
