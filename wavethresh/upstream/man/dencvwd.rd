\name{dencvwd}
\alias{dencvwd}
\title{
Calculate variances of wavlet coefficients of a p.d.f.
}
\usage{
dencvwd(hrproj, filter.number=hrproj$filter$filter.number,
    family=hrproj$filter$family, type="wavelet", bc="zero",
    firstk=hrproj$klim, RetFather=TRUE, verbose=FALSE)
}
\arguments{
\item{hrproj}{Output from \code{\link{denproj}} with \code{covar=T} argument.}
\item{filter.number}{The filter number of the wavelet basis to be used.
	This argument should not be altered from the default, as it is tied
	to the \code{hrproj} argument}
\item{family}{The family of wavelets to use.
	This argument should not be altered.}
\item{type}{The type of decomposition to be performed. This argument should not be altered.}
\item{bc}{The type of boundary conditions to be used.
	For density estimation this should always be zero.}
\item{firstk}{The bounds on the translation index of the empirical scaling function coefficients.}
\item{RetFather}{Ignore this.}
\item{verbose}{If TRUE the function will be chatty. Note that comments are only availble for part of the algorithm, so might not be very enlightening.}
}
\description{
Calculates the variances of the empirical wavelet coefficients by performing a 2D wavelet decomposition on the covariance matrix of the empirical scaling function coefficients of the probability density function.
}
\details{
This function is basically \code{\link{imwd}} adapted to handle zero
boundary conditions, except that only the variances are returned,
i.e. the diagonals of the covariance matrices produced.
Note that this code is not very efficient.
The full covariance matrices of all levels of coefficients are calculated,
and then the diagonals are extracted.
}
\value{
An object of class \code{\link{wd.object}}, but the contents are not a
standard wavelet transform, ie the object is used to hold other information
which organisationally is arranged like a wavelet tranform, ie variances of
coefficients. 
}
\seealso{\code{\link{denproj}},\code{\link{imwd}}}
\examples{
# Simulate data from the claw density, find the 
# empirical scaling function coefficients and covariances and then decompose
# both to give wavelet coefficients and their variances.

data <- rclaw(100)
datahr <- denproj(data, J=8, filter.number=2,family="DaubExPhase", covar=TRUE)
data.wd <- denwd(datahr)
\dontrun{plotdenwd(data.wd, top.level=(datahr$res$J-1))}

datavar <- dencvwd(datahr)
\dontrun{plotdenwd(datavar, top.level=(datahr$res$J-1))} 
}
\author{David Herrick}
\keyword{smooth}
