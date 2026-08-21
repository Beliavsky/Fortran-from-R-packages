\name{none}
\alias{none}
\title{
Locfit Evaluation Structure
}
\usage{
none()
}
\description{
  \code{none()} is an evaluation structure for \code{\link{locfit.raw}()},
  specifying no evaluation points. Only the initial parametric fit is
  computed - this is the easiest and most efficient way to coerce
  Locfit into producing a parametric regression fit.
}

\examples{
data(ethanol, package="locfit")
# fit a fourth degree polynomial using locfit
fit <- locfit(NOx~E,data=ethanol,deg=4,ev=none())
plot(fit,get.data=TRUE)
}

\keyword{smooth}
