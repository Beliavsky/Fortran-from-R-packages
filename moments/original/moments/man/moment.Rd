\name{moment}
\alias{moment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Statistical Moments }
\description{
    This function computes the sample moment of specified order. 
}
\usage{
moment(x, order = 1, central = FALSE, absolute = FALSE, na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{x}{a numeric vector of data.}
 \item{order}{order of the moment to be computed}
 \item{central}{a logical value - if central moments are
   to be computed.}
 \item{absolute}{a logical value - if absolute moments are
   to be computed.}
 \item{na.rm}{a logical value - remove NA values?}
} 
\author{ Lukasz Komsta }
\examples{
set.seed(1234)
x <- rnorm(10)
moment(x)
moment(x,order=3,absolute=TRUE)
}
\keyword{univar}
