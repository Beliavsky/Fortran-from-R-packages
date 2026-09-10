\name{dof}
\alias{dof}
\title{Compute number of non-zero coefficients in wd object}
\description{
Compute number of non-zero coefficients in \code{\link{wd}} object
}
\usage{
dof(wd)
}
\arguments{
\item{wd}{A \code{wavelet decomposition} object (such as that returned by the \code{\link{wd}} function).}
}
\details{
Very simple function that counts the number of non-zero coefficients in a \code{\link{wd}} class object. 
}
\value{
An integer that represents the number of non-zero coefficients in the input \code{\link{wd}} object. 
}
\section{RELEASE}{Version 3.0 Copyright Guy Nason 1994 }
\seealso{
\code{\link{wd}}, \code{\link{wd.object}}, \code{\link{threshold}}, \code{\link{threshold.wd}}. 
}
\examples{
#
# Let's generate some purely random numbers!!
#
myrandom <- rnorm(512)
#
# Take the discrete wavelet transform
#
myrandomWD <- wd(myrandom)
#
# How many coefficients are non-zero?
#
dof(myrandomWD)
# [1] 512
#
# All of them were nonzero!
#
# Threshold it
#
myrandomWDT <- threshold(myrandomWD, policy="universal")
#
# Now lets see how many are nonzero
#
dof(myrandomWDT)
# [1] 8
#
# Wow so 504 of the coefficients were set to zero! Spooky!
#
}
\keyword{models}
\author{G P Nason}
