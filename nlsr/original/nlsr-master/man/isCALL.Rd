\name{isCALL}
\alias{isCALL}
\title{isCALL
 
Test if argument is a particular call}
\usage{
isCALL(x, name)
}
\arguments{
\item{x}{object to be tested}

\item{name}{The function to test for.}
}
\description{
Test if \code{x} is a call to the function given by \code{name}.

Used in \code{\link{newSimplification}} definitions.
}
\examples{
x <- quote(mean(1:10))
isCALL(x, "mean")
isCALL(x, "sd")
}
\keyword{internal}