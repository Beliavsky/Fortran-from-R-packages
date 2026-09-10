\name{mwd}
\alias{mwd}
\title{Discrete multiple wavelet transform (decomposition).}
\description{
This function performs the discrete multiple wavelet transform (DMWT). Using an adaption of Mallat's pyramidal algorithm. The DMWT gives vector wavelet coefficients. 
}
\usage{
mwd(data, prefilter.type = "default", filter.type = "Geronimo",
        bc ="periodic", verbose = FALSE)
}
\arguments{
\item{data}{A vector containing the data you wish to decompose. The length of this vector must be a power of 2 times the dimension of the DMWT (multiplicity of wavelets).}
\item{prefilter.type}{This chooses the method of preprocessing required. The arguments will depend on filter.type, but "default" will always work.}
\item{filter.type}{Specifies which multi wavelet filter to use, The options are "\code{Geronimo}" (dimension 2) or "\code{Donovan3}" (dimension 3). The latter has not been tested fully and may contain bugs. See the function \code{\link{mfilter.select}} for further details.} 
\item{bc}{specifies the boundary handling. If \code{bc=="periodic"} the default, then the function you decompose is assumed to be periodic on its interval of definition, if \code{bc=="symmetric"} then the function beyond its boundaries is assumed to be a symmetric reflection of the function in the boundary.}
\item{verbose}{Controls the printing of "informative" messages whilst the computations progress. Such messages are generally annoying so it is turned off by default.}
}
\details{
The code implements Mallat's pyramid algorithm adapted for multiple wavelets using Xia, Geronimo, Hardin and Suter, 1996. The method takes a data vector of length \code{2^J*M}, and preprocesses it. This has two effects, firstly it puts the data into matrix form and then filters it so that the DMWT can operate more efficiently Most of the technical details are similar to the single wavelet transform except for the matrix algebra considerations, and the prefiltering process. See Downie and Silverman (1998) for further details and how this transform can be used in a statistical context.
}
\value{
An object of class \code{\link{mwd}}. 
}
\author{
Tim Downie
}
\section{RELEASE}{Version 3.9.6 (Although Copyright Tim Downie 1996)}
\seealso{
\code{\link{accessC.mwd}}, \code{\link{accessD.mwd}}, \code{\link{draw.mwd}}, \code{\link{mfirst.last}}, \code{\link{mfilter.select}}, \code{\link{mwd.object}}, \code{\link{mwr}}, \code{\link{plot.mwd}}, \code{\link{print.mwd}}, \code{\link{putC.mwd}}, \code{\link{putD.mwd}}, \code{\link{summary.mwd}}, \code{\link{threshold.mwd}}, \code{\link{wd}}, \code{\link{wr.mwd}}. 
}
\examples{
#
# Generate some test data
#
test.data <- example.1()$y
\dontrun{ts.plot(test.data)}
#
# Decompose test.data with multiple wavelet transform and
# plot the wavelet coefficients
#
tdmwd <- mwd(test.data)
\dontrun{plot(tdmwd)}
#[1] 1.851894 1.851894 1.851894 1.851894 1.851894 1.851894 1.851894
#
# You should see a plot with wavelet coefficients like in
#\code{\link{plot.wd}} but at each coefficient position
# there are two coefficients in two different colours one for each of
# the wavelets at that position.
#
# Note the scale for each level is returned by the function.
}
\keyword{math}
