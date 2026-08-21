\name{locfit.quasi}
\alias{locfit.quasi}
\title{
Local Quasi-Likelihood with global reweighting.
}
\usage{
locfit.quasi(x, y, weights, ..., iter=3, var=abs)
}
\description{
  \code{locfit.quasi} assumes a specified mean-variance relation,
  and performs iterartive reweighted local regression under this
  assumption. This is appropriate for local quasi-likelihood models,
  and is an alternative to specifying a family such as \code{"qpoisson"}.

  \code{locfit.quasi} is designed as a front end
  to \code{\link{locfit.raw}} with data vectors, or as an intemediary
  between \code{\link{locfit}} and \code{\link{locfit.raw}} with a
  model formula. If you can stand the syntax, the second calling
  sequence above will be slightly more efficient than the third.
}
\arguments{
\item{x}{
  Either a \code{\link{locfit}} model formula or a numeric vector
  of the predictor variable.
}
\item{y}{
  If \code{x} is numeric, \code{y} gives the response variable.
}
\item{weights}{Case weights to use in the fitting.}
\item{...}{
  Other arguments to \code{\link{locfit.raw}}
}
\item{iter}{Number of EM iterations to perform}
\item{var}{
  Function specifying the assumed relation between the mean and variance.
}
}
\value{
\code{"locfit"} object.
}
\seealso{
\code{\link{locfit}},
\code{\link{locfit.raw}}
}
\keyword{smooth}
