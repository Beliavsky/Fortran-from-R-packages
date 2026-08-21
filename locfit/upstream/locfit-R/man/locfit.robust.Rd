\name{locfit.robust}
\alias{locfit.robust}
\title{
Robust Local Regression
}
\usage{
locfit.robust(x, y, weights, ..., iter=3)
}
\description{
  \code{locfit.robust} implements a robust local regression where
  outliers are iteratively identified and downweighted, similarly
  to the lowess method (Cleveland, 1979). The iterations and scale
  estimation are performed on a global basis.

  The scale estimate is 6 times the median absolute residual, while
  the robust downweighting uses the bisquare function. These are
  performed in the S code so easily changed.

  This can be interpreted as an extension of M estimation to local
  regression. An alternative extension (implemented in locfit via
  \code{family="qrgauss"}) performs the iteration and scale estimation
  on a local basis.
}
\arguments{
\item{x}{
  Either a \code{\link{locfit}} model formula or a numeric vector
  of the predictor variable.
}
\item{y}{
  If \code{x} is numeric, \code{y} gives the response variable.
}
\item{weights}{weights to use in the fitting.}
\item{...}{Other arguments to \code{\link{locfit.raw}}.}
\item{iter}{Number of iterations to perform}
}
\value{
\code{"locfit"} object.
}
\seealso{
\code{\link{locfit}},
\code{\link{locfit.raw}}
}
\references{
Cleveland, W. S. (1979).
  Robust locally weighted regression and smoothing scatterplots.
  J. Amer. Statist. Assn. 74, 829-836.
}
\keyword{smooth}
% Converted by Sd2Rd version 0.2-a5.
