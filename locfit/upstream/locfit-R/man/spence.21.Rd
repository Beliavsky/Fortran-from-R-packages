\name{spence.21}
\alias{spence.21}
\title{
Spencer's 21 point graduation rule.
}
\usage{
spence.21(y)
}
\description{
  Spencer's 21 point rule is a weighted moving average operation
  for a sequence of observations equally spaced in time. The
  average at time t depends on the observations at times t-11,...,t+11.

  Except for boundary effects, the function will reproduce
  polynomials up to degree 3.
}
\arguments{
  \item{y}{Data vector of observations at equally spaced points.}
}
\value{
  A vector with the same length as the input vector, representing
  the graduated (smoothed) values.
}
\examples{
 data(spencer)
yy <- spence.21(spencer$mortality)
plot(spencer$age, spencer$mortality)
lines(spencer$age, yy)
}

\seealso{
  \code{\link{spence.15}},
  \code{\link{spencer}},
}

\references{
  Spencer, J. (1904).
  On the graduation of rates of sickness and mortality.
  Journal of the Institute of Actuaries 38, 334-343.
}

\keyword{smooth}
% Converted by Sd2Rd version 0.2-a5.
