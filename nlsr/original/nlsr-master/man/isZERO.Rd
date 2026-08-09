\name{isZERO}
\alias{isZERO}
\alias{isONE}
\alias{isMINUSONE}
\title{Test for constants}
\usage{
isZERO(x)
isONE(x)
isMINUSONE(x)
}
\arguments{
\item{x}{object to be tested}
}
\description{
Test for the values \code{0}, \code{1} or \code{-1}.

Used in \code{\link{newSimplification}} definitions.
}
\value{
Returns \code{TRUE} if the argument is the appropriate scalar value.
}
\examples{
isZERO(0)
x <- quote(0*1)
isZERO(x) # This is `*`(0, 1), not a value
isZERO(eval(x))
x <- quote(-1)
isMINUSONE(x) # This is `-`(1), not a value
isMINUSONE(eval(x))
}
\keyword{internal}