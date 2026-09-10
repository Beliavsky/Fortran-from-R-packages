\name{rmget}
\alias{rmget}
\title{Search for existing ipndacw matrices. }
\description{
Returns the integer corresponding to the smallest order \code{\link{ipndacw}} matrix of greater than or equal to order than the order, J requested. 

Not really intended for user use. 
}
\usage{
rmget(requestJ, filter.number, family)
}
\arguments{
\item{requestJ}{A positive integer representing the order of the \code{\link{ipndacw}} matrix that is \emph{required}.}
\item{filter.number}{The index number of the wavelet used to build the \code{\link{ipndacw}} matrix that is required.}
\item{family}{The wavelet family used to build the \code{\link{ipndacw}} matrix that is required.} 
}
\details{
Some of the matrices computed by \code{\link{ipndacw}} take a long time to compute. Hence it is a good idea to store them and reuse them. 

This function is asked to find an \code{\link{ipndacw}} matrix of a particular order, \emph{filter.number} and \emph{family}. The function steps through all of the directories in the \code{search()} list collecting names of all \code{\link{ipndacw}} matrices having the same \emph{filter.number} and \emph{family} characteristics. It then keeps any names where the \emph{order} is larger than, or equal to, the requested order. This means that a suitable \code{\link{ipndacw}} matrix of the same or larger order is visible in one of the \code{search()} directories. The matrix name with the smallest \code{order} is selected and the \emph{order} of the matrix is returned. The routine that called this function can then \code{get()} the matrix and either use it "as is" or extract the top-left hand corner of it if \code{requestJ} is less than the order returned by this function. 

If no such matrix, as described by the previous paragraph, exists then this function returns \code{NULL}. 

This function calls the subsidiary routine \code{\link{firstdot}}. 
}
\value{
If a matrix of order larger than or equal to the requested order exists somewhere on the search path \emph{and} the \code{filter.number} and \code{\link{family}} is as specified then its order is returned. If more than one such matrix exists then the order of the smallest one larger than or equal to the requested one is returned.
 
If no such matrix exists the function returns NULL. 
}
\references{Nason, G.P., von Sachs, R. and Kroisandt, G. (1998). Wavelet processes and adaptive estimation of the evolutionary wavelet spectrum. \emph{Technical Report}, Department of Mathematics University of Bristol/ Fachbereich Mathematik, Kaiserslautern.} 
\section{RELEASE}{Version 3.9 Copyright Guy Nason 1998 }
\seealso{
\code{\link{firstdot}}, \code{\link{ipndacw}}, \code{\link{rmname}}. 

}
\examples{
#
# Suppose there are no matrices in the search path.
#
# Let's look for the matrix rm.4.1.DaubExPhase (Haar wavelet matrix of
# order 4)
#
rmget(requestJ=4, filter.number=1, family="DaubExPhase")
#NULL
#
# I.e. a NULL return code. So there were no suitable matrices.
#
#If we create two Haar ipndacw matrix of order 7 and 8
#
ipndacw(-7, filter.number=1, family="DaubExPhase")
ipndacw(-8, filter.number=1, family="DaubExPhase")
#
# Now let's repeat the earlier search
#
rmget(requestJ=4, filter.number=1, family="DaubExPhase")
#[1] 7
#
# So, as we the smallest Haar ipndacw matrix available larger than
# the requested order of 4 is "7".
#
}
\author{G P Nason}
\keyword{manip}
