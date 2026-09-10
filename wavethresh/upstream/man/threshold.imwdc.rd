\name{threshold.imwdc}
\alias{threshold.imwdc}
\title{Threshold two-dimensional compressed wavelet decomposition object}
\description{
This function provides various ways to threshold a \code{imwdc} class object. 
}
\usage{
\method{threshold}{imwdc}(imwdc, verbose=FALSE, ...)
}
\arguments{
\item{imwdc}{The two-dimensional compressed wavelet decomposition object that you wish to threshold.}
\item{verbose}{if TRUE then the function prints out informative messages as it progresses.}
\item{\dots}{other arguments passed to the \code{\link{threshold.imwd}} function to control the thresholding characteristics such as policy, type of thresholding etc.}
}
\details{
This function performs exactly the same function as \code{\link{threshold.imwd}} except is accepts objects of class \code{imwdc} rather than imwd. Indeed, this function physically calls the \code{\link{threshold.imwd}} function after using the \code{\link{uncompress}} function to convert the input \code{imwdc} object into a \code{\link{imwd}} object. 
}
\value{
An object of class \code{imwdc} if the compression option is supplied and set to TRUE, otherwise a \code{\link{imwd}} object is returned. In either case the returned object contains the thresholded coefficients. Note that if the \code{return.threshold} option is set to TRUE then the threshold values will be returned rather than the thresholded object. }

\references{The FDR code segments were kindly donated by Felix Abramovich. }
\section{RELEASE}{Version 3.6 Copyright Guy Nason and others 1997 }
\seealso{
\code{\link{imwd}}, \code{\link{imwd.object}}, \code{\link{imwdc.object}}, \code{\link{threshold}}, \code{\link{uncompress}}. 
}
\examples{
#
# See examples in \code{\link{threshold.imwd}}.
#
}
\keyword{nonlinear}
\keyword{smooth}
\author{G P Nason}
