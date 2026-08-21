\name{aicplot}
\alias{aicplot}
\title{
Compute an AIC plot.
}
\description{
  The \code{aicplot} function loops through calls to the \code{\link{aic}}
function (and hence to \code{\link{locfit}}), using a different
smoothing parameter for each call.
The returned structure contains the AIC statistic for each fit, and can
be used to produce an AIC plot.
}
\usage{
aicplot(..., alpha)
}
\arguments{
  \item{...}{ arguments to the \code{\link{aic}}, \code{\link{locfit}}
     functions.}
  \item{alpha}{ Matrix of smoothing parameters. The \code{aicplot} function
     loops through calls to \code{\link{aic}}, using each row of
     \code{alpha} as the smoothing parameter in turn. If \code{alpha}
     is provided as a vector, it will be converted to a one-column
     matrix, thus interpreting each component as a nearest neighbor
    smoothing parameter.}
}

\value{
  An object with class \code{"gcvplot"}, containing the smoothing
  parameters and AIC scores. The actual plot is produced using
  \code{\link{plot.gcvplot}}.
}

\examples{
data(morths)
plot(aicplot(deaths~age,weights=n,data=morths,family="binomial",
  alpha=seq(0.2,1.0,by=0.05)))
}

\seealso{
  \code{\link{locfit}},
  \code{\link{locfit.raw}},
  \code{\link{gcv}},
  \code{\link{aic}},
  \code{\link{plot.gcvplot}}
}

\keyword{htest}
