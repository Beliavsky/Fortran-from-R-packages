\name{wst.object}
\alias{wst.object}
\title{(Packet ordered) Nondecimated wavelet transform decomposition objects.}
\description{
These are objects of class \code{wst} They represent a decomposition of a function with respect to a set of (all possible) shifted wavelets. 
}
\value{
The following components must be included in a legitimate `wst' object. 
\item{wp}{a matrix containing the packet ordered non-decimated wavelet coefficients. Each row of the matrix contains coefficients with respect to a particular resolution level. There are \code{nlevelsWT(wst)+1} rows in the matrix. Row \code{nlevels(wst)+1} (the ``bottom'') row contains the ``original'' data used to produce the wavelet packet coefficients. Rows \code{nlevels(wst)} to row 1 contain coefficients at resolution levels \code{nlevels(wst)-1} to 0 (so the first row contains coefficients at resolution level 0).
 
The columns contain the coefficients with respect to packets. A different packet length exists at each resolution level. The packet length at resolution level \code{i} is given by \code{2^i}. However, the \code{\link{getpacket.wst}} function should be used to access individual packets from a \code{\link{wst}} object.}
\item{Carray}{A matrix of the same dimensions and format as \code{wp} but containing the father wavelet coefficients.}
\item{nlevelsWT}{The number of levels in the decomposition. If you raise 2 to the power of \code{nlevels} you get the number of data points used in the decomposition.}
\item{filter}{a list containing the details of the filter that did the decomposition (equivalent to the return value from the \code{\link{filter.select}} function).}
\item{date}{The date that the transform was performed or the wst was modified.} 
}
\details{
To retain your sanity we recommend that the coefficients from a \code{wst} object be extracted in one of two ways: 

\itemize{
\item{use \code{\link{getpacket.wst}} to obtain individual packets of either father or mother wavelet coefficients.}
\item{use \code{\link{accessD.wst}} to obtain all mother coefficients at a particular resolution level.}
\item{use \code{\link{accessC.wst}} to obtain all father coefficients at a particular resolution level.}
}
 
You can obtain the coefficients directly from the \code{wst$wp} component (mother) or \code{wst$Carray} component (father) but you have to understand their organization described above. 
}

\section{GENERATION}{
This class of objects is returned from the \code{\link{wst}} function which computes the \emph{packets-ordered} non-decimated wavelet transform (effectively all possible shifts of the standard discrete wavelet transform). 

Many other functions return an object of class \code{wst}. 
}

\section{METHODS}{
The wst class of objects has methods for the following generic functions: \code{\link{AvBasis}}, \code{\link{InvBasis}}, \code{\link{LocalSpec}}, \code{\link{MaNoVe}}, \code{\link{accessC}}, \code{\link{accessD}}, \code{\link{convert}}, \code{\link{draw}}. \code{\link{getpacket}}. \code{\link{image}}. \code{\link{nlevelsWT}}, \code{\link{nullevels}}, \code{\link{plot}}, \code{\link{print}}, \code{\link{putC}}, \code{\link{putD}}, \code{\link{putpacket}}, \code{\link{summary}}, \code{\link{threshold}}. 
}
\section{RELEASE}{
Version 3.5.3 Copyright Guy Nason 1994 
}
\seealso{
\code{\link{wst}}
}
\keyword{classes}
\keyword{smooth}
\keyword{nonlinear}
\author{G P Nason}
