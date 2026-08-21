\name{cpplot}
\alias{cpplot}
\title{
Compute a Cp plot.
}
\usage{
cpplot(..., alpha, sig2)
}
\description{
  The \code{cpplot} function loops through calls to the \code{\link{cp}}
function (and hence to \code{link{locfit}}), using a different
smoothing parameter for each call.
The returned structure contains the Cp statistic for each fit, and can
be used to produce an AIC plot.
}
\arguments{
  \item{...}{ arguments to the \code{\link{cp}}, \code{\link{locfit}}
     functions.}
  \item{alpha}{ Matrix of smoothing parameters. The \code{cpplot} function
     loops through calls to \code{\link{cp}}, using each row of
     \code{alpha} as the smoothing parameter in turn. If \code{alpha}
     is provided as a vector, it will be converted to a one-column
     matrix, thus interpreting each component as a nearest neighbor
    smoothing parameter.}
  \item{sig2}{ Residual variance. If not specified, the residual variance
    is computed using the fitted model with the fewest residual degrees
    of freedom.}
}

\value{
  An object with class \code{"gcvplot"}, containing the smoothing
  parameters and CP scores. The actual plot is produced using
  \code{\link{plot.gcvplot}}.
}

\examples{
data(ethanol)
plot(cpplot(NOx~E,data=ethanol,alpha=seq(0.2,1.0,by=0.05)))
}

\seealso{
  \code{\link{locfit}},
  \code{\link{locfit.raw}},
  \code{\link{gcv}},
  \code{\link{aic}},
  \code{\link{plot.gcvplot}}
}

\keyword{htest}
% Converted by Sd2Rd version 0.2-a5.
