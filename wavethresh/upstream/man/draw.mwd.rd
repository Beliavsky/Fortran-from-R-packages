\name{draw.mwd}
\alias{draw.mwd}
\title{Draws a wavelet or scaling function used to compute an `mwd' object }
\description{
Draws picture of one wavelet or scaling function associated with the multiple wavelet decomposition object. \code{\link{mwd.object}}. 
}
\usage{
\method{draw}{mwd}(mwd, phi = 0, psi = 0, return.funct = FALSE, \dots)
}
\arguments{
\item{mwd}{The \code{\link{mwd}} class object whose associated wavelet or scaling function you wish to draw. }
\item{phi}{description not yet available}
\item{psi}{If \code{phi} is non-zero then the `phi'-th scaling function of the wavelet family used for mwd will be plotted. \code{phi} must be between 0 and \code{mwd$filter$nphi}.}
\item{return.funct}{If true then the vector used as phi/psi in the plot command is returned.}
\item{\dots}{Additional arguments to pass to the \code{plot} function}
}
\details{
It is usual to specify just one of phi and psi. IF neither phi nor psi are specified then phi=1 is the default. An error is generated if both phi=0 and psi=0 or if both are nonzero. 
}
\value{
If the \code{return.funct} argument is set to \code{TRUE} then the function values in the plot are returned otherwise \code{NULL} is returned. 
}
\note{If the \code{return.funct} argument is \code{FALSE} a plot of the mother wavelet or scaling function is plotted on the active graphics device. 
}
\section{RELEASE}{Version 3.9.6 (Although Copyright Tim Downie 1995-6).}
\seealso{
\code{\link{accessC.mwd}}, \code{\link{accessD.mwd}}, \code{\link{mfirst.last}}, \code{\link{mfilter.select}}, \code{\link{mwd}}, \code{\link{mwd.object}}, \code{\link{mwr}}, \code{\link{plot.mwd}}, \code{\link{print.mwd}}, \code{\link{putC.mwd}}, \code{\link{putD.mwd}}, \code{\link{summary.mwd}}, \code{\link{threshold.mwd}}, \code{\link{wd}}, \code{\link{wr.mwd}}. 
}
\examples{
#
# Do a multiple wavelet decomposition on vector: ynoise
#
ynoise <- rnorm(512, sd = 0.1)
ymwd <- mwd(ynoise,filter.type="Geronimo")
#
# Draw a picture of the second Geronimo wavelet.
#
\dontrun{draw(ymwd,psi=2)}
#
#
}
\keyword{hplot}
\author{G P Nason}
