\name{irregwd.objects}
\alias{irregwd.objects}
\title{Irregular wavelet decomposition objects.}
\description{
These are objects of classes 

\code{wd}

They represent a decomposition of a function with respect to a wavelet basis. The function will have been interpolated to a grid and these objects represent the discrete wavelet transform \code{\link{wd}}.
}
\section{GENERATION}{
This class of objects is returned from the \code{\link{irregwd}} function. Some other functions that process these kinds of objects also return this class of object (such as \code{\link{threshold.irregwd}}.)
}
\section{METHODS}{
The \code{irregwd} class of objects has methods for the following generic functions: \code{\link{plot}}, \code{\link{threshold}},
}
\section{STRUCTURE}{
All components in a legitimate `irregwd' are identical to the components in an ordinary \code{\link{wd.object}} with the exception of \code{type} component and with the addition of the following component: 
\describe{
\item{c}{vector that aids in the calculation of variances of wavelet coefficients (used by \code{\link{threshold.irregwd}}).}
}
}
\details{
To retain your sanity the C and D coefficients should be extracted by the \code{\link{accessC}} and \code{\link{accessD}} functions and inserted using the \code{\link{putC}} and \code{\link{putD}} functions (or more likely, their methods), rather than by the \code{$} operator.
 
One can use the \code{\link{accessc}} function to obtain the \code{c} component. 

Mind you, if you want to muck about with coefficients directly, then you'll have to do it yourself by working out what the fl.dbase list means (see \code{\link{first.last}} for a description.) 
}
\section{RELEASE}{
Version 3.9.4 Copyright Arne Kovac 1997, Help Copyright Guy Nason 2004 
}
\seealso{
\code{\link{irregwd}}, \code{\link{threshold.irregwd}}, \code{\link{plot.irregwd}},\code{\link{wd}} 
}
\keyword{smooth} 
\author{G P Nason}
