\name{AutoBasis}
\alias{AutoBasis}
\title{Run Coifman-Wickerhauser best basis algorithm on wavelet packet object}
\description{
Runs the Coifman-Wickerhauser best basis algorithm on a wavelet packet
object. Packets not in the basis are replaced by vectors of NAs.
Superceded by the \code{\link{MaNoVe}} functions.
}
\details{
Superceded by the \code{\link{MaNoVe}} functions (which run in C code).
}
\value{A wp class object which contains the select basis. All packets that are
	not in the basis get replaced by vectors of NAs.
}
\seealso{\code{\link{MaNoVe}}}
\author{G P Nason}
\keyword{math}
