\name{linfnorm}
\alias{linfnorm}
\title{Compute L infinity distance between two vectors of numbers. }
\description{
Compute L infinity distance between two vectors of numbers (maximum absolute difference between two vectors). 
}
\usage{
linfnorm(u,v)
}
\arguments{
\item{u}{first vector of numbers} 
\item{v}{second vector of numbers} 
}
\details{
Function simply computes the L infinity distance between two vectors and is implemented as 

\code{max(abs(u-v))}
}
\value{
A real number which is the L infinity distance between two vectors. 
}
\note{
This function would probably be more accurate if it used the Splus function \code{vecnorm}.} 
\section{RELEASE}{Version 3.6 Copyright Guy Nason 1995 }
\seealso{
\code{\link{l2norm}}, \code{\link{wstCV}}, \code{\link{wstCVl}}. 
}
\examples{
#
# What is the L infinity norm between the following sets of vectors
#
p <- c(1,2,3,4,5)
q <- c(1,2,3,4,5)
r <- c(2,3,4,5,6)
linfnorm(p,q)
# [1] 0
linfnorm(q,r)
# [1] 1
linfnorm(r,p)
# [1] 1
}
\keyword{algebra}
\author{G P Nason}
