\name{catch}
\alias{catch}
\title{Catching Results, Warnings and Errors Simultaneously}
\description{
   Catches results, warnings and errors.
}
\usage{
catch(expr)
}
\arguments{
  \item{expr}{expression to be evaluated, typically a function call.}
}
\value{
  \code{\link{list}} with components:
  \item{\code{value}}{value of \code{expr} or \code{NULL} in case of
    an error.}
  \item{\code{warning}}{warning message (see \code{\link{simpleWarning}} or
    \code{\link{warning}()}) or \code{\link{NULL}} in case of no warning.}
  \item{\code{error}}{error message (see \code{\link{simpleError}} or
    \code{\link{stop}()}) or \code{\link{NULL}} in case of no error.}
}
\details{
  This function is particularly useful for large(r) simulation studies
  to not fail until finished.
}
\author{
  Marius Hofert (based on \code{doCallWE()} and \code{tryCatch.W.E()} in
  the \R package \pkg{simsalapar}).
}
\examples{
catch(log(2))
catch(log(-1))
catch(log("a"))
}
\keyword{programming}
