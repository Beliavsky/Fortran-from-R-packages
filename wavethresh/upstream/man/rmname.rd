\name{rmname}
\alias{rmname}
\title{Return a ipndacw matrix style name.}
\description{
This function returns a character string according to a particular format for naming \code{\link{ipndacw}} matrices. 
}
\usage{
rmname(J, filter.number, family)
}
\arguments{
\item{J}{A negative integer representing the order of the \code{\link{ipndacw}} matrix.} 
\item{filter.number}{The index number of the wavelet used to build the \code{\link{ipndacw}} matrix.} 
\item{family}{The wavelet family used to build the \code{\link{ipndacw}} matrix.} 
}
\details{
Some of the matrices computed by \code{\link{ipndacw}} take a long time to compute. Hence it is a good idea to store them and reuse them. This function generates a name according to a particular naming scheme that permits a search algorithm to easily find the matrices. 

Each matrix has three defining characteristics: its \emph{order}, \emph{filter.number} and \emph{family}. Each of these three characteristics are concatenated together to form a name. 
}
\value{
A character string containing the name of a matrix according to a particular naming scheme. 
}
\references{Nason, G.P., von Sachs, R. and Kroisandt, G. (1998). Wavelet processes and adaptive estimation of the evolutionary wavelet spectrum. \emph{Technical Report}, Department of Mathematics University of Bristol/ Fachbereich Mathematik, Kaiserslautern. }
\section{RELEASE}{Version 3.9 Copyright Guy Nason 1998 }
\seealso{
\code{\link{ewspec}}, \code{\link{ipndacw}}, 
}
\examples{
#
# What's the name of the order 4 Haar matrix?
#
rmname(-4, filter.number=1, family="DaubExPhase")
#[1] "rm.4.1.DaubExPhase"
#
# What's the name of the order 12 Daubechies least-asymmetric wavelet
# with 7 vanishing moments?
#
rmname(-12, filter.number=7, family="DaubLeAsymm")             
#[1] "rm.12.7.DaubLeAsymm"
}
\keyword{manip}
\author{G P Nason}
