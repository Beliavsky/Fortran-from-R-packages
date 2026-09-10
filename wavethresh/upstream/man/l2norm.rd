\name{l2norm}
\alias{l2norm}
\title{Compute L2 distance between two vectors of numbers. }
\description{
Compute L2 distance between two vectors of numbers (square root of sum of squares of differences between two vectors). 
}
\usage{
l2norm(u,v)
}
\arguments{
\item{u}{first vector of numbers} 
\item{v}{second vector of numbers} 
}
\details{
Function simply computes the L2 distance between two vectors and is implemented as 

\code{sqrt(sum((u-v)^2))}
}
\value{
A real number which is the L2 distance between two vectors. 
}
\note{This function would probably be more accurate if it used the Splus function \code{vecnorm}.} 
\section{RELEASE}{Version 3.6 Copyright Guy Nason 1995}
\seealso{
\code{\link{linfnorm}}, \code{\link{wstCV}}, \code{\link{wstCVl}}. 
}
\examples{
#
# What is the L2 norm between the following sets of vectors
#
p <- c(1,2,3,4,5)
q <- c(1,2,3,4,5)
r <- c(2,3,4,5,6)
l2norm(p,q)
# [1] 0
l2norm(q,r)
# [1] 2.236068
l2norm(r,p)
# [1] 2.236068
}
\keyword{algebra}
\author{G P Nason}
