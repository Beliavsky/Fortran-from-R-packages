\name{lscvplot}
\alias{lscvplot}
\title{
Compute the LSCV plot.
}
\usage{
lscvplot(..., alpha)
}
\description{
  The \code{lscvplot} function loops through calls to the \code{\link{lscv}}
function (and hence to \code{link{locfit}}), using a different
smoothing parameter for each call.
The returned structure contains the LSCV statistic for each density
estimate, and can be used to produce an LSCV plot.
}
\arguments{
  \item{...}{ arguments to the \code{\link{lscv}}, \code{\link{locfit}}
     functions.}
  \item{alpha}{ Matrix of smoothing parameters. The \code{lscvplot} function
     loops through calls to \code{\link{lscv}}, using each row of
     \code{alpha} as the smoothing parameter in turn. If \code{alpha}
     is provided as a vector, it will be converted to a one-column
     matrix, thus interpreting each component as a nearest neighbor
    smoothing parameter.}
}

\value{
  An object with class \code{"gcvplot"}, containing the smoothing
  parameters and LSCV scores. The actual plot is produced using
  \code{\link{plot.gcvplot}}.
}

\seealso{
  \code{\link{locfit}},
  \code{\link{locfit.raw}},
  \code{\link{gcv}},
  \code{\link{lscv}},
  \code{\link{plot.gcvplot}}
}

\keyword{htest}
% Converted by Sd2Rd version 0.2-a5.
