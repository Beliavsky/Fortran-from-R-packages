\name{print.nv}
\alias{print.nv}
\title{Print a node vector object, also used by several other functions to
	obtain packet list information}
\usage{
\method{print}{nv}(x, printing = TRUE, verbose = FALSE, ...)
}
\arguments{
\item{x}{The \code{\link{nv.object}} that you wish to print}
\item{printing}{If FALSE then nothing is printed. This argument is
	here because the results of the printing are also useful
	to many other routines where you want the results but are not
	bothered by actually seeing the results}
\item{verbose}{Not actually used}
\item{\dots}{Other arguments}
}
\description{
Ostensibly prints out node vector information, but also produces packet
	indexing information for several functions.
}
\details{
A node vector contains selected basis information, but this is stored
as a tree object. Hence, it is not immediately obvious which basis elements
have been stored. This function produces a list of the packets at each
resolution level that have been selected in the basis. This information
is so useful to other functions that the function is used even when
printing is not the primary objective.
}
\value{
A list containing two components: \code{indexlist} and \code{rvector}.
The former is a list of packets that were selected at each resolution
level. Rvector encodes a list of "rotate/non-rotate" instructions in binary.
At each selected packet level a decision has to be made whether to select
the LH or RH basis element, and this information is stored in \code{rvector}.
}
\seealso{
\code{\link{InvBasis.wst}},
\code{\link{nv.object}},
\code{\link{plot.wp}}}
\examples{
v <- rnorm(128)
vwst <- wst(v)
vnv <- MaNoVe(vwst)
print(vnv)
#Level :  6  Action is  R (getpacket Index:  1 )
#Level :  5  Action is  L (getpacket Index:  2 )
#Level :  4  Action is  L (getpacket Index:  4 )
#Level :  3  Action is  R (getpacket Index:  9 )
#Level :  2  Action is  L (getpacket Index:  18 )
#There are  6  reconstruction steps
#
# The L or R indicate whether to move to the left or the right basis function
# when descending the node tree
#
#
}
\author{G P Nason}
\keyword{print}
