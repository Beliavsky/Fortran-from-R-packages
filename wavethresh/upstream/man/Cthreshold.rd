\name{Cthreshold}
\alias{Cthreshold}
\title{Calls C code to threshold wd class object.}
\usage{
Cthreshold(wd, thresh.type = "soft", value = 0, levels = 3:(nlevelsWT(wd) - 1)) 
}
\arguments{
\item{wd}{The wavelet object that you wish to threshold.}
\item{thresh.type}{The type of thresholding. This can be "soft" or "hard". See \code{\link{threshold}} and methods for further details.}
\item{value}{The threshold value that you want to be used (e.g. for hard thresholding wavelet coefficients whose absolute value is less than}
\item{levels}{The resolution levels that you wish to compute the threshold on and apply the threshold to.}
}
\description{
A routine that calls a C code function to do thresholding. This is really a test routine to call a C thresholding function (Cthreshold) and the user is advised to use the R based generic thresholding function
\code{\link{threshold}} and/or its methods as they contain a wider range of thresholding options.
}
\details{
For general use it is recommended to use the
\code{\link{threshold}} functions as they have a wider variety of options
and also work for more complex varieties of wavelet transforms
(i.e. non-decimated, complex-valued, etc).
However, in the right, limited, situation this function can be useful.
This function directly calls the C thresholding function Cthreshold().
The C function is used by routines that operate on behalf of the function
that carries out two-fold cross validation in C (\code{\link{CWCV}}) which
is also accessible using the \code{policy="cv"} option too
\code{\link{threshold.wd}}

This function can be used by the user. It might be a bit faster than
\code{\link{threshold.wd}} but mostly because it is simpler and does less checking than
\code{\link{threshold.wd}}.
}
\value{
A \code{\link{wd.object}} class object, but containing thresholded coefficients.
}
\seealso{\code{\link{threshold}}}
\examples{
#
# See copious examples in the help to threshold.wd
#
}
\author{G P Nason}
\keyword{smooth}
