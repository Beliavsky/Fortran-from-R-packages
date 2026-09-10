\name{denwd}
\alias{denwd}
\title{Wavelet decomposition of empirical scaling function coefficients of a p.d.f.  }
\usage{
denwd(coef)
}
\arguments{
\item{coef}{Output from \code{\link{denproj}}}
}
\description{
Performs wavelet decomposition on the empirical scaling function coefficients of the probability density function.
}
\details{
The empirical scaling function coefficients are decomposed using the DWT with zero boundary conditions.
}
\value{
An object of class \code{\link{wd.object}}
}
\seealso{\code{\link{denproj}},\code{\link{plotdenwd}},\code{\link{wd}},
	\code{\link{denwr}}}
\examples{
# Simulate data from the claw density, find the empirical
# scaling function coefficients, decompose them and plot
# the resulting wavelet coefficients
 
data <- rclaw(100)
datahr <- denproj(data, J=8, filter.number=2,family="DaubExPhase")
data.wd <- denwd(datahr)
\dontrun{plotdenwd(data.wd, top.level=(datahr$res$J-1))}
}
\author{David Herrick}
\keyword{smooth}
