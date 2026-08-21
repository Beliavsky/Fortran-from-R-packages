\name{lcvplot}
\alias{lcvplot}
\title{
Compute the likelihood cross-validation plot.
}
\usage{
lcvplot(..., alpha)
}
\description{
  The \code{lcvplot} function loops through calls to the \code{\link{lcv}}
function (and hence to \code{link{locfit}}), using a different
smoothing parameter for each call.
The returned structure contains the likelihood cross validation statistic
for each fit, and can be used to produce an LCV plot.
}
\arguments{
  \item{...}{ arguments to the \code{\link{lcv}}, \code{\link{locfit}}
     functions.}
  \item{alpha}{ Matrix of smoothing parameters. The \code{aicplot} function
     loops through calls to \code{\link{lcv}}, using each row of
     \code{alpha} as the smoothing parameter in turn. If \code{alpha}
     is provided as a vector, it will be converted to a one-column
     matrix, thus interpreting each component as a nearest neighbor
    smoothing parameter.}
}

\value{
  An object with class \code{"gcvplot"}, containing the smoothing
  parameters and LCV scores. The actual plot is produced using
  \code{\link{plot.gcvplot}}.
}

\examples{
data(ethanol)
plot(lcvplot(NOx~E,data=ethanol,alpha=seq(0.2,1.0,by=0.05)))
}

\seealso{
  \code{\link{locfit}},
  \code{\link{locfit.raw}},
  \code{\link{gcv}},
  \code{\link{lcv}},
  \code{\link{plot.gcvplot}}
}

\keyword{htest}
% Converted by Sd2Rd version 0.2-a5.
