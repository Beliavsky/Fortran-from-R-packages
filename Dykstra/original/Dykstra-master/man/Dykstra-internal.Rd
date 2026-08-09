\name{Dykstra-internal}
\alias{print.dykstra}
\title{
Print Dykstra Results
}
\description{
Prints convergence and fit information.
}
\details{
Displays the convergence flag, the final iteration count, the obtained solution, and the value of the qudratic function at the obtained solution.
}

\usage{
\method{print}{dykstra}(x, n = 5L, ...)
}
\arguments{
  \item{x}{
  	Object of class "dykstra" output from \code{\link{dykstra}}.
}
  \item{n}{
  	Does not print solutions with more than \code{n} coefficients.
}
  \item{\ldots}{
  Ignored.
}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\examples{
?dykstra
}
\keyword{ internal }