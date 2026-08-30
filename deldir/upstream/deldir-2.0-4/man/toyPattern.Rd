\name{toyPattern}
\alias{toyPattern}
\docType{data}
\title{
    A toy marked point pattern object, with 59 points.
}
\description{
    A simulated object of class \code{"ppp"} provided for use
    in an example illustrating the application of \code{deldir()}
    to point pattern objects. 
}
\usage{toyPattern}
\format{
  An object of class \code{"ppp"} consisting of a simulated marked
  point pattern.  Entries include
  \tabular{ll}{
    \code{x} \tab Cartesian \eqn{x}-coordinates \cr
    \code{y} \tab Cartesian \eqn{y}-coordinates \cr
    \code{marks} \tab factor with levels \code{"a","b","c","d"}
    }
}
\source{ Simulated.  }
\examples{
dtoy    <- deldir(toyPattern) # "Tags" are the marks of the pattern.
set.seed(42)
dtoy.nt <- deldir(toyPattern,z=round(runif(59),2)) # Tags are numeric.
}
\keyword{datasets}
