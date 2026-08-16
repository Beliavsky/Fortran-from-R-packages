\name{predict.polynomial}
\alias{predict.polynomial}
\title{Evaluate a Polynomial}
\usage{
\method{predict}{polynomial}(object, newdata, \dots)
}
\arguments{
  \item{object}{A polynomial object to be evaluated.}
  \item{newdata}{Argument at which evaluation is requested.  May be
    numeric or itself a polynomial}
  \item{\dots}{Not used by this method.}
}
\description{
  Evaluate a polynomial at a given numeric or polynomial argument.
}
\details{
  This is a method for the generic function \code{\link{predict}}.

  The polynomial is evaluated according to the Horner scheme for speed
  and numerical accuracy.
}
\value{
  Evaluated object of the same class as \code{newdata}.
}
\seealso{
  \code{\link{as.function.polynomial}}
}
\keyword{symbolmath}
