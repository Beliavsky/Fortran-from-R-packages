\name{tpwd}
\alias{tpwd}
\title{Tensor product 2D wavelet transform}
\usage{
tpwd(image, filter.number = 10, family = "DaubLeAsymm", verbose = FALSE)
}
\arguments{
\item{image}{The image  you wish to subject to the tensor product WT}
\item{filter.number}{The smoothness of wavelet, see \code{\link{filter.select}}}
\item{family}{The wavelet family you wish to use}
\item{verbose}{Whether or not you wish to print out informative messages}
}
\description{
Performs the tensor product 2D wavelet transform. This is a
related, but different, 2D wavelet transform compared to
\code{\link{imwd}}.
}
\details{
The transform works by first taking the regular 1D wavelet
transform across all columns in the image and storing these coefficients
line by line back into the image. Then to this new image we apply
the regular 1D wavelet transform across all rows in the image.

Hence, the top-left coefficient is the smoothed version both horizontally
and vertically. The left-most row contains the image smoothed horiztonally,
but then detail picked up amongst the horizontal smooths vertically.

Suggested by Rainer von Sachs.
}
\value{
A list with the following components:
\item{tpwd}{A matrix with the same dimensions as the input \code{image},
	but containing the tensor product wavelet transform coefficients.}
\item{filter.number}{The filter number used}
\item{family}{The wavelet family used}
\item{type}{The type of transform used}
\item{bc}{The boundary conditions used}
\item{date}{When the transform occurred}
}
\seealso{\code{\link{imwd}},\code{\link{tpwr}}}
\examples{
data(lennon)
ltpwd <- tpwd(lennon)
\dontrun{image(log(abs(ltpwd$tpwd)), col=grey(seq(from=0, to=1, length=100)))}
}
\author{G P Nason}
\keyword{math}
