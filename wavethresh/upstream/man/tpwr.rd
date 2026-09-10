\name{tpwr}
\alias{tpwr}
\title{Inverse tensor product 2D wavelet transform.}
\usage{
tpwr(tpwdobj, verbose = FALSE)
}
\arguments{
\item{tpwdobj}{An object which is a list which contains the items
	indicated in the return value of \code{\link{tpwd}}}
\item{verbose}{Whether informative messages are printed}
}

\description{
Performs the inverse transform to \code{\link{tpwd}}.
}
\details{
Performs the inverse transform to \code{\link{tpwd}}.
}
\value{
A matrix, or image, containing the inverse tensor product wavelet
transform of the image contained in the \code{tpwd} component of the
\code{tpwdobj} object.
}
\seealso{\code{\link{imwr}},\code{\link{tpwd}}}
\examples{
data(lennon)
ltpwd <- tpwd(lennon)
#
# now perform the inverse and compare to the original
#
ltpwr <- tpwr(ltpwd)
sum((ltpwr - lennon)^2)
# [1] 9.22802e-10
}
\author{G P Nason}
\keyword{math}
