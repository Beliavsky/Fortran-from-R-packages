\name{neldermead.destroy}
\alias{neldermead.destroy}
\title{
  Erase a neldermead object.
}
\description{
  \code{neldermead.destroy} calls \code{optimbase.destroy} and
  \code{optimsimplex.destroy} to erase the content of \code{this$optbase} and
  \code{this$simplex0}.
}
\usage{
  neldermead.destroy(this = NULL)
}
\arguments{
  \item{this}{A neldermead object.}
}
\value{
  Return an updated neldermead object.
}
\seealso{
  \code{\link[optimbase]{optimbase.destroy}},
  \code{\link[optimsimplex]{optimsimplex.destroy}}
}
\author{
Author of Scilab neldermead module: Michael Baudin (INRIA - Digiteo)

Author of R adaptation: Sebastien Bihorel (\email{sb.pmlab@gmail.com})
}
\keyword{ method }
