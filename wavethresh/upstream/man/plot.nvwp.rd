\name{plot.nvwp}
\alias{plot.nvwp}
\title{Depict wavelet packet basis specfication}
\usage{
\method{plot}{nvwp}(x, \dots)
}
\arguments{
\item{x}{The wavelet packet node vector you wish to plot, nvwp class object}
\item{\dots}{Other arguments to the central plot function} 
}
\description{
The nvwp class object (generated from \code{\link{MaNoVe.wp}} for example)
contains a wavelet packet basis specification. This function produces
a graphical depiction of such a basis.
}
\details{
The vertical axis indicates the resolution level, the horizontal axes
indicates the packet index for the finest scales.
}
\value{
Nothing
}
\seealso{\code{\link{MaNoVe.wp}},\code{\link{print.nvwp}},\code{\link{wp}}}
\examples{
v <- rnorm(512)
vwp <- wp(v)
vnv <- MaNoVe(vwp)
\dontrun{plot(vnv)}
}
\author{G P Nason}
\keyword{hplot}
