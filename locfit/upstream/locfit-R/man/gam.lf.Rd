\name{gam.lf}
\alias{gam.lf}
\title{
Locfit call for Generalized Additive Models
}
\usage{
gam.lf(x, y, w, xeval, ...)
}
\description{
  This is a locfit calling function used by
  \code{\link{lf}()} terms in additive models. It is
  not normally called directly by users.
}
\arguments{
  \item{x}{numeric predictor}
  \item{y}{numeric response}
  \item{w}{prior weights}
  \item{xeval}{evaluation points}
  \item{...}{other arguments to \code{\link{locfit.raw}()}}
}

\seealso{
  \code{\link{locfit}},
  \code{\link{locfit.raw}},
  \code{\link{lf}},
  \code{gam}
}

\keyword{models}
