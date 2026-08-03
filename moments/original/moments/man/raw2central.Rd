\name{raw2central}
\alias{raw2central}
\title{ Raw to central moments }
\description{
  This function transforms a vector, matrix or data frame of raw moments to
  a vector, matrix or data frame of central moments.
}
\usage{
raw2central(mu.raw)
}
\arguments{
  \item{mu.raw}{ A numeric vector, matrix or data frame of raw moments.  For a vector, mu.raw[0] is
                 the order 0 raw moment, mu.raw[1] is the order 1 raw moment and so forth.
                 For a matrix or data frame, row vector mu.raw[0,] contains the order 0 raw moments,
                 row vector mu.raw[1,] contains the order 1 raw moments and so forth. }
}
\value{
  A vector matrix or data frame of central moments.  For matrices and data frame, column vectors
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
  \code{\link{central2raw}}
}
\examples{
set.seed(1234)
x <- rnorm(10000)
mu.raw.x <- all.moments( x, order.max=4 )
mu.central.x <- all.moments( x, central=TRUE, order.max=4 )
raw2central( mu.raw.x )
mu.central.x
M <- matrix( x, nrow=1000, ncol=10 )
mu.raw.M <- all.moments( M, order.max=4 )
mu.central.M <- all.moments( M, central=TRUE, order.max=4 )
raw2central( mu.raw.M )
mu.central.M
D <- data.frame( M )
mu.raw.D <- all.moments( D, order.max=4 )
mu.central.D <- all.moments( D, central=TRUE, order.max=4 )
raw2central( mu.raw.D )
mu.central.D
}
\keyword{ math }
\keyword{ univar }
