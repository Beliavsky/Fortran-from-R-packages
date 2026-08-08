\name{neldermead.search}
\alias{neldermead.search}
\title{
  Starts the optimization
}
\description{
  Performs the optimization associated with the method associated with the
  \code{method} element of the neldermead object and find the optimum. If the
  \code{restartflag} element is enabled, automatic restarts are performed, based
  on the \code{restartdetection} element.
}
\usage{
  neldermead.search(this = NULL)
}
\arguments{
  \item{this}{A neldermead object.}
}
\value{
  Return an updated neldermead object.
}
\author{
Author of Scilab neldermead module: Michael Baudin (INRIA - Digiteo)

Author of R adaptation: Sebastien Bihorel (\email{sb.pmlab@gmail.com})
}
\seealso{
  \code{\link[neldermead]{fminsearch}},
  \code{\link[neldermead]{neldermead}},
  \code{\link[neldermead]{neldermead.set}},
}
\keyword{ method }
