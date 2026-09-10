\name{wp.object}
\alias{wp.object}
\title{Wavelet Packet decomposition objects.}
\description{
These are objects of classes \code{wp} They represent a decomposition of a function with respect to a set of wavelet packet functions. 
}
\details{
To retain your sanity we recommend that wavelet packets be extracted in one of two ways: 

\itemize{
\item{use \code{\link{getpacket.wp}} to obtain individual packets.}
\item{use \code{\link{accessD.wp}} to obtain all coefficients at a particular resolution level.}
}
 
You can obtain the coefficients directly from the \code{wp$wp} component but you have to understand their organization described above. 
}
\value{
The following components must be included in a legitimate `wp' object. 
\item{wp}{a matrix containing the wavelet packet coefficients. Each row of the matrix contains coefficients with respect to a particular resolution level. There are \code{nlevelsWT(wp)+1} rows in the matrix. Row \code{nlevels(wp)+1} (the ``bottom'') row contains the ``original'' data used to produce the wavelet packet coefficients. Rows \code{nlevels}(wp) to row 1 contain coefficients at resolution levels \code{nlevels(wp)-1} to 0 (so the first row contains coefficients at resolution level 0). 

The columns contain the coefficients with respect to packets. A different packet length exists at each resolution level. The packet length at resolution level \code{i} is given by \code{2^i}. However, the \code{\link{getpacket.wp}} function should be used to access individual packets from a \code{\link{wp}} object.}
\item{nlevelsWT}{The number of levels in the wavelet packet decomposition. If you raise 2 to the power of nlevels you get the number of data points used in the decomposition. }
\item{filter}{a list containing the details of the filter that did the decomposition (equivalent to the return value from the \code{\link{filter.select}} function). }
\item{date}{The date that the transform was performed or the wp was modified.} 
}
\section{GENERATION}{
This class of objects is returned from the \code{\link{wp}} function to represent a wavelet packet decomposition of a function. Many other functions return an object of class wp. 
}
\section{METHODS}{
The wp class of objects has methods for the following generic functions: \code{\link{InvBasis}}, \code{\link{MaNoVe}}, \code{\link{accessC}}, \code{\link{accessD}}, \code{\link{basisplot}}, \code{\link{draw}}. \code{\link{getpacket}}, \code{\link{nlevelsWT}}, \code{\link{plot}}, \code{\link{print}}, \code{\link{putC}}, \code{\link{putD}}, \code{\link{putpacket}}, \code{\link{summary}}, \code{\link{threshold}}. 
}
\section{RELEASE}{
Version 3.5.3 Copyright Guy Nason 1994 
}
\seealso{
\code{\link{wp}}
}
\keyword{classes}
\keyword{smooth}
\author{G P Nason}
