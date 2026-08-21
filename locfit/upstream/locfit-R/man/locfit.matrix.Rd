\name{locfit.matrix}
\alias{locfit.matrix}
\title{
Reconstruct a Locfit model matrix.
}
\usage{
locfit.matrix(fit, data)
}
\description{
  Reconstructs the model matrix, and associated variables such as
  the response, prior weights and censoring indicators, from a
  \code{locfit} object. This is used by functions such as
  \code{\link{fitted.locfit}}; it is not normally called directly.
  The function will only work properly if the data frame has not been
  changed since the fit was constructed.
}
\arguments{
\item{fit}{Locfit object}
\item{data}{ Data Frame.}
}
%\item{...}{
%Other arguments to \code{\link{locfit.raw}()} (or the \code{lfproc}).
%}
\value{
A list with variables \code{x} (the model matrix); \code{y} (the response);
\code{w} (prior weights); \code{sc} (scales); \code{ce} (censoring indicator)
and \code{base} (baseline fit).
}
\seealso{
\code{\link{locfit}}, \code{\link{fitted.locfit}}, \code{\link{residuals.locfit}}
}
\keyword{models}
