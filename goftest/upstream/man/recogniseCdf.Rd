\name{recogniseCdf}
\alias{recogniseCdf}
\title{
  Explanatory Name of Distribution Function
}
\description{
  Recognises many standard cumulative distribution functions and
  returns a string describing the distribution.
}
\usage{
recogniseCdf(s="punif")
}
\arguments{
  \item{s}{
    A single character string giving the \emph{name}
    of an \R function that calculates cumulative
    probabilities.
  }
}
\details{
  The list of recognised distribution functions
  includes all those available in the \pkg{stats} package
  and in \pkg{goftest}.

  By convention, the name of a cumulative distribution function
  begins with the letter \code{p}. For example, \code{punif} is the
  cumulative distribution function of the uniform distribution.

  The initial letter \code{p} can be omitted in the
  function \code{recogniseCdf}.
}
\value{
  Character string, or \code{NULL} if the name is not recognised.
}
\author{
  Adrian Baddeley.
}
\seealso{
  \code{\link{pAD}}
}
\examples{
   recogniseCdf("punif")
   recogniseCdf("unif")
   recogniseCdf("pt")
}
\keyword{distribution}
