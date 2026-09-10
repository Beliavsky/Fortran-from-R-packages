\name{Psiname}
\alias{Psiname}
\title{Return a PsiJ list object style name.}
\description{
This function returns a character string according to a particular format for naming \code{\link{PsiJ}} objects. 
}
\usage{
Psiname(J, filter.number, family)
}
\arguments{
\item{J}{A negative integer representing the order of the \code{\link{PsiJ}} object.}
\item{filter.number}{The index number of the wavelet used to build the \code{\link{PsiJ}} object.}
\item{family}{The wavelet family used to build the \code{\link{PsiJ}} object.}
}
\details{
Some of the objects computed by \code{\link{PsiJ}} take a long time to compute. Hence it is a good idea to store them and reuse them. This function generates a name according to a particular naming scheme that permits a search algorithm to easily find the matrices. 

Each object has three defining characteristics: its \emph{order}, \emph{filter.number} and \emph{family}. Each of these three characteristics are concatenated together to form a name. 

This function performs exactly the same role as \code{\link{rmname}} except for objects produced by \code{\link{PsiJ}}. 
}
\value{
A character string containing the name of an object according to a particular naming scheme. 
}
\references{Nason, G.P., von Sachs, R. and Kroisandt, G. (1998). Wavelet processes and adaptive estimation of the evolutionary wavelet spectrum. \emph{Technical Report}, Department of Mathematics University of Bristol/ Fachbereich Mathematik, Kaiserslautern}
\section{RELEASE}{Version 3.9 Copyright Guy Nason 1998 }
\seealso{
\code{\link{PsiJ}}
}
\examples{
#
# What's the name of the order 4 Haar PsiJ object?
#
Psiname(-4, filter.number=1, family="DaubExPhase")
#[1] "Psi.4.1.DaubExPhase"
#
# What's the name of the order 12 Daubechies least-asymmetric wavelet PsiJ
# with 7 vanishing moments?
#
Psiname(-12, filter.number=7, family="DaubLeAsymm")             
#[1] "Psi.12.7.DaubLeAsymm"
}
\keyword{character}
\author{G P Nason}
