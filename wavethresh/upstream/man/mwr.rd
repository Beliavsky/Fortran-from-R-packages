\name{mwr}
\alias{mwr}
\title{Multiple discrete wavelet transform (reconstruction). }
\description{
This function performs the reconstruction stage of Mallat's pyramid algorithm adapted for multiple wavelets (see Xia et al.(1996)), i.e. the discrete inverse \emph{multiple} wavelet transform. 
}
\usage{
mwr(mwd, prefilter.type = mwd$prefilter, verbose = FALSE, start.level = 0,
	returnC = FALSE)
}
\arguments{
\item{mwd}{A multiple wavelet decomposition object as returned by \code{\link{mwd}}.} 
\item{prefilter.type}{Usually best not to change this (i.e. not to use a different prefilter on the reconstruction to the one used on decomposition).}
\item{verbose}{Controls the printing of "informative" messages whilst the computations progress. Such messages are generally annoying so it is turned off by default.}
\item{start.level}{The level you wish to start reconstruction at. The is usually the first (level 0).} 
\item{returnC}{If this is FALSE then a vector of the same length as the argument data supplied to the function \code{\link{mwd}} that constructed the supplied \code{\link{mwd.object}}. is returned, Ie. the reconstructed data. If true then the last level (highest resolution) C coefficients are returned in matrix form. This matrix has not been postprocessed. }
}
\details{
The code implements Mallat's pyramid algorithm adapted for multiple wavelet decompositions (Xia et al. 1996). In the reconstruction the quadrature mirror filters G and H are supplied with C0 and D0, D1, ... D(J-1) (the wavelet coefficients) and rebuild C1,..., CJ. 

The matrix CJ is postprocessed which returns the full reconstruction 

If \code{\link{mwd.object}} was obtained directly from \code{\link{mwd}} then the original function can be reconstructued exactly. Usually, the \code{\link{mwd.object}} has been modified in some way, for examples, some coefficients set to zero by \code{\link{threshold}}. Mwr then reconstructs the function with that set of wavelet coefficients. 

See also Downie and Silverman, 1998 
}
\value{
Either a vector containing the final reconstruction or a matrix containing unpostprocessed coefficients. 
}
\section{RELEASE}{Version 3.9.6 (Although Copyright Tim Downie 1996)}
\seealso{
\code{\link{accessC.mwd}}, \code{\link{accessD.mwd}}, \code{\link{draw.mwd}}, \code{\link{mfirst.last}}, \code{\link{mfilter.select}}, \code{\link{mwd}}, \code{\link{mwd.object}}, \code{\link{plot.mwd}}, \code{\link{print.mwd}}, \code{\link{putC.mwd}}, \code{\link{putD.mwd}}, \code{\link{summary.mwd}}, \code{\link{threshold.mwd}}, \code{\link{wd}}, \code{\link{wr.mwd}}. 
}
\examples{
#
# Decompose and then exactly reconstruct test.data
#
test.data <- rnorm(128)
tdecomp <- mwd(test.data)
trecons <- mwr(tdecomp)
#
# Look at accuracy of reconstruction
max(abs(trecons - test.data))
#[1] 2.266631e-12
#
# See also the examples of using \code{\link{wr}} or mwr in
# the \code{examples} section of
# the help for \code{\link{threshold.mwd}}.
}
\keyword{manip}
\author{Tim Downie}
