\name{firstdot}
\alias{firstdot}
\title{Return the location of the first period character within a character string (for a vector of strings of arbitrary length). }
\description{
Returns the index of the location of the first period character within a character string for a series of strings in a vector of character string of arbitrary length). 

This is a subsidiary routine for \code{\link{rmget}} and not really intended for user use. 
}
\usage{
firstdot(s)
}
\arguments{
\item{s}{Vector of character strings.}
}
\details{
A very simple function. It searches through a character string for the first period character and the returns the position of that period character. It performs this search for each of the character strings in the input vector. 
}
\value{
A vector of integers of the same length as the input vector. Each integer in the output vector is the index position of the first period character in the corresponding character string in the input vector. If a character string does not contain a period character then the corresponding output integer is zero. 
}
\references{Nason, G.P., von Sachs, R. and Kroisandt, G. (1998). Wavelet processes and adaptive estimation of the evolutionary wavelet spectrum. \emph{Technical Report}, Department of Mathematics University of Bristol/ Fachbereich Mathematik, Kaiserslautern.}
\section{RELEASE}{Version 3.9 Copyright Guy Nason 1998 }
\seealso{
\code{\link{rmget}}
}
\examples{
#
# Let's find the first dot in the following strings...
#
firstdot("mary.had.a.little.lamb")
#[1] 5
#
# I.e. the first period was after "mary" -- the fifth character
#
# This following string doesn't have any periods in it.
#
firstdot("StellaArtois")
#[1] 0
#
# The function works on vectors of character strings
#
TopCricketAve <- c("Don.Bradman", "Graeme.Pollock", "George.Headley",
	"Herbert.Sutcliffe", "Vinod.Kambli", "Javed.Miandad")
firstdot(TopCricketAve)
#[1] 4 7 7 8 6 6
}
\keyword{utilities}
\author{G P Nason}
