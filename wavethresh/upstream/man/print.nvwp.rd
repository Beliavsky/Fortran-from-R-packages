\name{print.nvwp}
\alias{print.nvwp}
\title{Print a wavelet packet node vector object,
	also used by several other functions to
	obtain packet list information}
\usage{
\method{print}{nvwp}(x, printing = TRUE, ...)
}
\arguments{
\item{x}{The nvwp that you wish to print}
\item{printing}{If FALSE then nothing is printed. This argument is
	here because the results of the printing are also useful
	to many other routines where you want the results but are not
	bothered by actually seeing the results}
\item{\dots}{Other arguments}
}
\description{
Ostensibly prints out wavlet packet node vector information,
	but also produces packet
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
A list containing two components: \code{level} and \code{pkt}.
These are the levels and packet indices of the select packets in the basis.
}
\seealso{
\code{\link{InvBasis.wp}},
\code{\link{MaNoVe.wp}},
\code{\link{plot.nvwp}},
\code{\link{plot.wp}}}
\examples{
v <- rnorm(128)
vwp <- wp(v)
vnv <- MaNoVe(vwp)
print(vnv)
#Level:  6  Packet:  1
#Level:  3  Packet:  0
#Level:  2  Packet:  4
#Level:  2  Packet:  13
#Level:  2  Packet:  15
#Level:  1  Packet:  5
#Level:  1  Packet:  10
#Level:  1  Packet:  13
#Level:  1  Packet:  14
#Level:  1  Packet:  15
#Level:  1  Packet:  16
#Level:  1  Packet:  20
#Level:  1  Packet:  21
#Level:  1  Packet:  24
#Level:  0  Packet:  8
#Level:  0  Packet:  9
#Level:  0  Packet:  12
#Level:  0  Packet:  13
#Level:  0  Packet:  14
#Level:  0  Packet:  15
#Level:  0  Packet:  22
#Level:  0  Packet:  23
#Level:  0  Packet:  24
#Level:  0  Packet:  25
#Level:  0  Packet:  34
#Level:  0  Packet:  35
#Level:  0  Packet:  36
#Level:  0  Packet:  37
#Level:  0  Packet:  38
#Level:  0  Packet:  39
#Level:  0  Packet:  44
#Level:  0  Packet:  45
#Level:  0  Packet:  46
#Level:  0  Packet:  47
#Level:  0  Packet:  50
#Level:  0  Packet:  51
#Level:  0  Packet:  56
#Level:  0  Packet:  57
#Level:  0  Packet:  58
#Level:  0  Packet:  59
}
\author{G P Nason}
\keyword{print}
