\name{c2to4}
\alias{c2to4}
\usage{
c2to4(index)
}
\arguments{
\item{index}{The integer you wish to convert}
}
\title{Take integer, represent in binary, then think of and return that
	representation in base 4}
\description{
Not designed, or really useful, for casual user use!

For example: take the integer 5. In binary this is 101. Then, this
representation in base 4 is 16+1 =17.

This function is used by \code{\link{accessD.wpst}} to help it access
coefficients.
}
\details{
Description says all
}
\value{
The converted number
}
\seealso{\code{\link{accessD.wpst}}}
\examples{
c2to4(5)
}
\author{G P Nason}
\keyword{misc}
