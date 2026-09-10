\name{nv.object}
\alias{nv.object}
\title{Node vector objects.}
\description{
These are objects of classes 

\code{nv}

They represent a basis in a packet-ordered non-decimated wavelet transform object. 
}
\details{
A \code{nv} object is a description of a basis which is a path through a packet ordered non-decimated wavelet transform. To view the basis just print it! See the examples in \code{\link{numtonv}} for a print out of its structure.
 
A similar object exists for describing a basis in a wavelet packet object see nvwp. 
}
\value{
The following components must be included in a legitimate `nv' object. 

\item{node.list}{This is a complicated structure composed of one-dimensional array of \code{nv$nlevelsWT} lists. Each item in the array is itself a list having two components\code{$upperctrl} and \code{upperl}. Each component is described as follows:
\describe{
\item{upperctrl}{The `upperctrl' item in each is the most important. It consists of a vector of characters. Each character refers to a node in the non-decimated wavelet tree at that level and can only be one of the characters L (for left), R (for right) and S (for stop). Each character in the vector informs reconstruction algorithms that, to do the best thing (whatever the best thing is in any particular case, e.g. select the minimum entropy node downwards), you should select the left/right node or stop at the current node.}
\item{upperl}{The `upperl' vector is in 1-1 correspondance with the `upperctrl' vector. Each entry is a number related in some way to the L/R/S entry. (For the minumum entropy this is the minmum entropy achieved by this selection).}

\item{nlevelsWT}{The number of levels in the \code{\link{wst}} object that was involved in the creation of the \code{nv} object. Nv objects describe a basis relative to a packet ordered non-decimated wavelet transform object and thus must know the number of levels in that object.}
}
}
}
\section{GENERATION}{
This class of objects is returned from the \code{\link{MaNoVe.wst}} and \code{\link{numtonv}} functions. The former returns the minimum entropy basis (most sparse basis) obtained using the Coifman-Wickerhauser, 1992 algorithm. The latter permits selection of a basis by an index number. 
}
\section{METHODS}{
The \code{nv} class of objects has methods for the following generic functions: print, \code{\link{nlevelsWT}}, \code{\link{InvBasis}}, 
}
\section{RELEASE}{
Version 3.6.0 Copyright Guy Nason 1995
}
\seealso{
\code{\link{wst}}, \code{\link{wst.object}}, \code{\link{numtonv}}, \code{\link{print}}, \code{\link{nlevelsWT}}, \code{\link{InvBasis}}, \code{\link{MaNoVe.wst}}. 
}
\keyword{classes}
\author{G P Nason}
