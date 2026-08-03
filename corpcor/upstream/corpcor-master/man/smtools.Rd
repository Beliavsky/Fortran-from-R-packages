\name{smtools}
\alias{sm2vec}
\alias{sm.index}
\alias{vec2sm}

\title{Some Tools for Handling Symmetric Matrices}

\usage{
sm2vec(m, diag = FALSE)
sm.index(m, diag = FALSE)
vec2sm(vec, diag = FALSE, order = NULL)
}

\arguments{
  \item{m}{symmetric matrix}
  \item{diag}{logical. Should the diagonal be included in the conversion to and from a vector?}
  \item{vec}{vector of unique elements from a symmetric matrix}
  \item{order}{order of the entries in \code{vec}}
}
\description{
  \code{sm2vec} takes a symmetric matrix and puts
  the lower triagonal entries into a vector (cf. \code{\link{lower.tri}}).

  \code{sm.index} lists the corresponding x-y-indices for each entry
  in the vector produced by \code{sm2vec}.
  
  \code{vec2sm} reverses the operation by  \code{sm2vec} and converts the
  vector back to a symmetric matrix. If \code{diag=FALSE} the
  diagonal of the resulting matrix will consist of NAs.  If \code{order}
  is supplied then the input vector \code{vec} will first be rearranged accordingly.
 }

\value{
  A vector (\code{sm2vec}), a two-column matrix with indices (\code{sm.index}),
  or a symmetric matrix (\code{vec2sm}).
}


\author{
  Korbinian Strimmer (\url{https://strimmerlab.github.io/}).
}

\seealso{\code{\link{lower.tri}}.}

\examples{
# load corpcor library
library("corpcor")

# a symmetric matrix
m = rbind(
 c(3,1,1,0),
 c(1,3,0,1),
 c(1,0,2,0),
 c(0,1,0,2)
)
m

# convert into vector (including the diagonals)
v = sm2vec(m, diag=TRUE)
v.idx = sm.index(m, diag=TRUE)
v
v.idx

# put back to symmetric matrix
vec2sm(v, diag=TRUE)

# convert from vector with specified order of the elements
sv = sort(v)
sv
ov = order(v)
ov
vec2sm(sv, diag=TRUE, order=ov)
}
\keyword{utilities}
