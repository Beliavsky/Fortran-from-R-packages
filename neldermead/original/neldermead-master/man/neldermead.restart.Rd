\name{neldermead.restart}
\alias{neldermead.restart}
\title{
  Restart neldermead search.
}
\description{
  Update the simplex with \code{neldermead.updatesimp} and restart the search
  with \code{neldermead.search}.
}
\usage{
  neldermead.restart(this = NULL)
}
\arguments{
  \item{this}{A neldermead object.}
}
\value{
  Returns an updated neldermead object.
}
\author{
Author of Scilab neldermead module: Michael Baudin (INRIA - Digiteo)

Author of R adaptation: Sebastien Bihorel (\email{sb.pmlab@gmail.com})
}
\seealso{
  \code{\link[neldermead]{neldermead.updatesimp}},
  \code{\link[neldermead]{neldermead.search}},
}
\keyword{ method }
