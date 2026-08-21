\name{cpar}
\alias{cpar}
\title{
Conditionally parametric term for a Locfit model.
}
\usage{
cpar(x,...)
}
\description{
  A term entered in a \code{\link{locfit}} model formula using
  \code{cpar} will result in a fit that is conditionally parametric.
  Equivalent to \code{lp(x,style="cpar")}.

  This function is presently almost deprecated. Specifying a conditionally
  parametric fit as \code{y~x1+cpar(x2)} wil no longer work; instead, the
  model is specified as \code{y~lp(x1,x2,style=c("n","cpar"))}.
}

\arguments{
  \item{x}{numeric variable.}
  \item{...}{Other arguments to \code{link{lp}()}.}
}

\examples{
data(ethanol, package="locfit")
# fit a conditionally parametric model
fit <- locfit(NOx ~ lp(E, C, style=c("n","cpar")), data=ethanol)
plot(fit)
# one way to force a parametric fit with locfit
fit <- locfit(NOx ~ cpar(E), data=ethanol)
}

\seealso{
  \code{\link{locfit}}
}

\keyword{models}
