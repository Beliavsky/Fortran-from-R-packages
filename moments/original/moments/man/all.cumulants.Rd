\name{all.cumulants}
\alias{all.cumulants}
\title{ Statistical Cumulants }
\description{
  This function calculates the cumulants for all orders specified in 
  the given vector, matrix or data frame of raw moments}
\usage{
all.cumulants(mu.raw)
}
\arguments{
  \item{mu.raw}{ A numeric vector, matrix or data frame of raw moments.  For a vector, mu.raw[0] is
                 the order 0 raw moment, mu.raw[1] is the order 1 raw moment and so forth.
                 For a matrix or data frame, row vector mu.raw[0,] contains the order 0 raw moments,
                 row vector mu.raw[1,] contains the order 1 raw moments and so forth. }
}
\details{
  No details are currently described.
}
\value{
  A vector matrix or data frame of cumulants.  For matrices and data frame, column vectors
  correspond to different random variables.
}
\references{
Papoulis, A., Pillai, S. U. (2002) Probability, Random Variables and Stochastic Processes, Fourth Edition,
  McGraw-Hill, New York, 146-147.

}
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\seealso{
  \code{\link{moment}},
  \code{\link{all.moments}}
}
\examples{
set.seed(1234)
x <- rnorm(10000)
mu.raw.x <- all.moments( x, order.max=6 )
all.cumulants( mu.raw.x )
M <- matrix( x, nrow=1000, ncol=10 )
mu.raw.M <- all.moments( M, order.max=6 )
all.cumulants( mu.raw.M )
D <- data.frame( M )
mu.raw.D <- all.moments( D, order.max=6 )
all.cumulants( mu.raw.D )
}
\keyword{ math }
\keyword{ univar }
