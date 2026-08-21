\name{HotellingsTsq}
\alias{qHotellingsTsq}
\alias{pHotellingsTsq}
\title{Hotellings T square distribution}
\description{
  The Hotellings T square distribution is the distribution of the squared
  Mahalanobis distances with respected to estimated variance covariance
  matrices. 
}
\usage{
qHotellingsTsq(p,n,m)
pHotellingsTsq(q,n,m)
          }
\arguments{
  \item{p}{a (vector of) probabilities}
  \item{q}{a vector of quantils}
  \item{n}{number of parameters, the p parameter of Hotellings
    \eqn{T^2}{T^2} distribution}
  \item{m}{number of dimensions, the m parameter of the Hotellings \eqn{T^2}{T^2}
    distribution}
}
\value{
  \item{qHotellingsT2}{a vector of quantils}
  \item{pHotellingsT2}{a vector of probabilities}
}

\details{
  The Hotellings \eqn{T^2} with paramter p and m is the distribution
  empirical squared Mahalanobis distances of a m dimensional vector with respect
  to a  variance covariance matrix estimated based on np degrees of freedom. 
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{ellipses}}, \code{\link{ConfRadius}}, \code{\link{pf}} 
}
\examples{
(q <- qHotellingsTsq(seq(0,0.9,by=0.1),3,25))
pHotellingsTsq(q,3,25)
}
\keyword{multivariate}
