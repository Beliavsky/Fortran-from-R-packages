\name{lflim}
\alias{lflim}
\title{
Construct Limit Vectors for Locfit fits.
}
\usage{
lflim(limits, nm, ret)
}
\description{
  This function is used internally to interpret \code{xlim} and \code{flim}
  arguments. It should not be called directly.
}
\arguments{
\item{limits}{
Limit argument.
}
\item{nm}{
Variable names.
}
\item{ret}{
Initial return vector.
}
}
\value{
Vector with length 2*dim.
}
\seealso{
\code{\link{locfit}}
}
\keyword{smooth}
