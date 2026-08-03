\name{central2raw}
\alias{central2raw}
\title{ Central to raw moments }
\description{
  This function transforms a vector, matrix or data frame of central moments to
  a vector, matrix or data frame of raw moments.
}
\usage{
central2raw(mu.central,eta)
}
\arguments{
  \item{mu.central}{ A numeric vector, matrix or data frame of central moments.  For a vector, mu.central[0] is
                     the order 0 central moment, mu.central[1] is the order 1 central moment and so forth.
                     For a matrix or data frame, row vector mu.central[0,] contains the order 0 central moments,
                     row vector mu.central[1,] contains the order 1 central moments and so forth. }
  \item{eta}{ A numeric vector of sample mean or expected values }                     
}
\value{
  A vector matrix or data frame of raw moments.  For matrices and data frame, column vectors
  correspond to different random variables.
}
\references{
Papoulis, A., Pillai, S. U. (2002) Probability, Random Variables and Stochastic Processes, Fourth Edition,
  McGraw-Hill, New York, 146-147.

}
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\seealso{
  \code{\link{moment}},
  \code{\link{all.moments}},
  \code{\link{raw2central}}
}
\examples{
set.seed(1234)
x <- rnorm(10000)
mu.raw.x <- all.moments( x, order.max=4 )
eta.x <- mu.raw.x[2]
mu.central.x <- all.moments( x, central=TRUE, order.max=4 )
central2raw( mu.central.x, eta.x )
mu.raw.x
M <- matrix( x, nrow=1000, ncol=10 )
mu.raw.M <- all.moments( M, order.max=4 )
eta.M <- mu.raw.M[2,]
mu.central.M <- all.moments( M, central=TRUE, order.max=4 )
central2raw( mu.central.M, eta.M )
mu.raw.M
D <- data.frame( M )
mu.raw.D <- all.moments( D, order.max=4 )
eta.D <- mu.raw.D[2,]
mu.central.D <- all.moments( D, central=TRUE, order.max=4 )
central2raw( mu.central.D, eta.D )
mu.raw.D
}
\keyword{ math }
\keyword{ univar }
