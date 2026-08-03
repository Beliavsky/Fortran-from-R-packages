\name{all.moments}
\alias{all.moments}
\title{ Statistical Moments }
\description{
  This function computes all the sample moments of the chosen type up to a given order.
}
\usage{
all.moments(x, order.max = 2, central = FALSE, absolute = FALSE, na.rm = FALSE)
}
\arguments{
  \item{x}{ A numeric vector, matrix or data frame of data.  For matrices and data frames, each column is a random variable }
  \item{order.max}{ the maximum order of the moments to be computed with a default value of 2. }
  \item{central}{ a logical value, if TRUE, central moments are computed.  Otherwise, raw moments are computed }
  \item{absolute}{ a logical value, if TRUE, absolute moments are computed.  Otherwise, standard moments are computed }
  \item{na.rm}{ a logical value, if TRUE, remove NA values.  Otherwise, keep NA values }
}
\details{
  The minimum value for order.max is 2.  The function stops running for values less than 2 and
  the message "maximum order whould be at least 2" is displayed on standard output.
}
\value{
  A vector, matrix or data frame of moments depending on the nature of the argument x.  If x is a vector,
  then the value returned is a vector, say mu, where mu[1] is the order 0 moment, mu[2] is the order 1 moment
  and so forth.  If x is a matrix or data frame, then the value returned is a matrix or data frame, respectively.
  In this case, suppose mu is the value returned.  Then, row vector mu[1,] contains the order 0 moments,
  mu[2,] contains the order 1 moments and so forth.
}
\references{
Papoulis, A., Pillai, S. U. (2002) Probability, Random Variables and Stochastic Processes, Fourth Edition,
  McGraw-Hill, New York, 146-147.

}
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\seealso{
  \code{\link{moment}},
  \code{\link{raw2central}}
}
\examples{
set.seed(1234)
x <- rnorm(10000)
all.moments( x, order.max=4 )
all.moments( x, central=TRUE, order.max=4 )
all.moments( x, absolute=TRUE, order.max=4 )
all.moments( x, central=TRUE, absolute=TRUE, order.max=4 )
M <- matrix( x, nrow=1000, ncol=10 )
all.moments( M, order.max=4 )
all.moments( M, central=TRUE, order.max=4 )
all.moments( M, absolute=TRUE, order.max=4 )
all.moments( M, central=TRUE, absolute=TRUE, order.max=4 )
D <- data.frame( M )
all.moments( D, order.max=4 )
all.moments( D, central=TRUE, order.max=4 )
all.moments( D, absolute=TRUE, order.max=4 )
all.moments( D, central=TRUE, absolute=TRUE, order.max=4 )
}
\keyword{ math }
\keyword{ univar }
