\name{tryDebugger}
\alias{tryDebugger}
\title{Empirical variograms for compositions}
\description{
  An R-debugger that also works with errors in parameters.
}
\usage{
tryDebugger(dump = last.dump)
}
\arguments{
  \item{dump}{An R dump object created by 'dump.frames'.}
}
\value{
  Nothing. 
}
\details{
  Works like debugger, with the small exception that it also works in
  situations of nasty errors, like recursive parameter evaluation,
  missing parameters, and additional errors in arguments. 
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link[utils]{debugger}}
}

\examples{
\dontrun{
f <- function(x,y=y) {y}
f(1)
tryDebugger() # works
debugger() # Does not allow to browse anything
}
}
\keyword{debugging}
