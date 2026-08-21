\name{hatmatrix}
\alias{hatmatrix}
\title{
Weight diagrams and the hat matrix for a local regression model.
}
\usage{
hatmatrix(formula, dc=TRUE, \dots)
}
\arguments{
  \item{formula}{model formula.}
  \item{dc}{derivative adjustment (see \code{\link{locfit.raw}})}
  \item{...}{Other arguments to \code{\link{locfit}} and
\code{\link{locfit.raw}}.
}}
\description{
\code{hatmatrix()} computes the weight diagrams (also known as
equivalent or effective kernels) for a local regression smooth.
Essentially, \code{hatmatrix()} is a front-end to \code{\link{locfit}()},
setting a flag to compute and return weight diagrams, rather than the
fit.
}
\value{
A matrix with n rows and p columns; each column being the
weight diagram for the corresponding \code{locfit} fit point.
If \code{ev="data"}, this is the transpose of the hat matrix.
}
\seealso{
\code{\link{locfit}}, \code{\link{plot.locfit.1d}}, \code{\link{plot.locfit.2d}},
\code{\link{plot.locfit.3d}}, \code{\link{lines.locfit}}, \code{\link{predict.locfit}}
}
\keyword{smooth}
