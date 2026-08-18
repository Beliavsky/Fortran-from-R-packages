\name{colMedians}
\alias{colMedians}
\title{Medians of Columns}

\description{
   Compute the sample medians of the columns (non-rows) of a data.frame or array.
}

\usage{
colMedians( x, na.rm = FALSE )
}

\arguments{
   \item{x}{a data.frame or array.}
   \item{na.rm}{a logical value indicating whether \code{NA} values should be
      stripped before the computation proceeds.}
}

\value{
   A vector or array of the medians of each column (non-row) of \code{x}
   with dimension \code{dim( x )[-1]}.
}

\seealso{\code{\link{rowMedians}},\code{\link{median}},\code{\link{colMeans}}.}

\author{Arne Henningsen}

\examples{
   data( "Electricity", package = "Ecdat" )
   colMedians( Electricity )

   a4 <- array( 1:120, dim = c(5,4,3,2),
      dimnames = list( c("a","b","c","d","e"), c("A","B","C","D"),
      c("x","y","z"), c("Y","Z") ) )
   colMedians( a4 )
   median( a4[ , "B", "x", "Z" ] )  # equal to
   colMedians( a4 )[ "B", "x", "Z" ]
}

\keyword{array}
