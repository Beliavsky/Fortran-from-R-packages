\name{wd.dh}
\alias{wd.dh}
\title{Compute specialized wavelet transform for density estimation}
\usage{
wd.dh(data, filter.number = 10, family = "DaubLeAsymm", type = "wavelet",
    bc = "periodic", firstk = NULL, verbose = FALSE)
}
\arguments{
\item{data}{The father wavelet coefficients}
\item{filter.number}{The smoothness of the underlying wavelet to use,	
	see \code{\link{filter.select}}}
\item{family}{The wavelet family to use, see \code{\link{filter.select}}}
\item{type}{The type of wavelet to use}
\item{bc}{Type of boundarie conditions}
\item{firstk}{A parameter that originates from \code{\link{denproj}}}
\item{verbose}{If \code{TRUE} then informative messages are printed.}
}
\description{
Computes the discrete wavelet transform, but with zero boundary conditions
especially for density estimation.
}
\details{
This is a subsidiary routine, not intended for direct user use for density
estimation. The main routines for wavelet density estimation are
\code{\link{denwd}}, \code{\link{denproj}}, \code{\link{denwr}}.

The input to this function should be projected father wavelet coefficients
as computed by \code{\link{denproj}}, but usually supplied to this function
by \code{\link{denwd}}.

Thresholding should be carried out by the user independently of these
functions.
}
\seealso{
\code{\link{denproj}},
\code{\link{denwd}}}
\value{
An object of class \code{\link{wd}}, but assumed on the basis of
zero boundary conditions.
}
\author{David Herrick}
\keyword{math}
\keyword{smooth}
