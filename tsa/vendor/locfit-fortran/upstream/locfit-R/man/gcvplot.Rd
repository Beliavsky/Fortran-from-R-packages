\name{gcvplot}
\alias{gcvplot}
\title{
Compute a generalized cross-validation plot.
}
\usage{
gcvplot(..., alpha, df=2)
}
\description{
  The \code{gcvplot} function loops through calls to the \code{\link{gcv}}
function (and hence to \code{link{locfit}}), using a different
smoothing parameter for each call.
The returned structure contains the GCV statistic for each fit, and can
be used to produce an GCV plot.
}
\arguments{
  \item{...}{ arguments to the \code{\link{gcv}}, \code{\link{locfit}}
     functions.}
  \item{alpha}{ Matrix of smoothing parameters. The \code{gcvplot} function
     loops through calls to \code{\link{gcv}}, using each row of
     \code{alpha} as the smoothing parameter in turn. If \code{alpha}
     is provided as a vector, it will be converted to a one-column
     matrix, thus interpreting each component as a nearest neighbor
    smoothing parameter.}
  \item{df}{ Degrees of freedom to use as the x-axis. 2=trace(L), 
    3=trace(L'L).}
}

\value{
  An object with class \code{"gcvplot"}, containing the smoothing
  parameters and GCV scores. The actual plot is produced using
  \code{\link{plot.gcvplot}}.
}

\examples{
data(ethanol)
plot(gcvplot(NOx~E,data=ethanol,alpha=seq(0.2,1.0,by=0.05)))
}

\seealso{
  \code{\link{locfit}},
  \code{\link{locfit.raw}},
  \code{\link{gcv}},
  \code{\link{plot.gcvplot}},
  \code{\link{summary.gcvplot}}
}

\keyword{htest}
% Converted by Sd2Rd version 0.2-a5.
