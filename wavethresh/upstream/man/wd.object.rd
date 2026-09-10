\name{wd.object}
\alias{wd.object}
\title{Wavelet decomposition objects}
\description{
These are objects of classes 

\code{wd}

They represent a decomposition of a function with respect to a wavelet basis (or tight frame in the case of the (time-ordered) non-decimated wavelet decomposition). 
}
\value{
The following components must be included in a legitimate `wd' object. 

\item{C}{a vector containing each level's smoothed data. The wavelet transform works by applying both a smoothing filter and a bandpass filter to the previous level's smoothed data. The top level contains data at the highest resolution level. Each of these levels are stored one after the other in this vector. The matrix 
	
\code{fl.dbase$first.last.c}
	
determines exactly where each level is stored in the vector. Likewise, coefficients stored when the NDWT has been used should only be extracted using the ``access'' and ``put'' functions below.}
\item{D}{wavelet coefficients. If you were to write down the discrete wavelet transform of a function then these D would be the coefficients of the wavelet basis functions. Like the C, they are also formed in a pyramidal manner, but stored in a linear array. The storage details are to be found in 

\code{fl.dbase$first.last.d}

Likewise, coefficients stored when the NDWT has been used should only be extracted using the ``access'' and ``put'' functions below. }
\item{nlevelsWT}{The number of levels in the pyramidal decomposition that produces the coefficients. If you raise 2 to the power of nlevels you get the number of data points used in the decomposition.}
\item{fl.dbase}{The first last database associated with this decomposition. This is a list consisting of 2 integers, and 2 matrices. The matrices detail how the coefficients are stored in the C and D components of the `wd.object'. See the help on \code{\link{first.last}} for more information. }
\item{filter}{a list containing the details of the filter that did the decomposition}
\item{type}{either \bold{wavelet} indicating that the ordinary wavelet transform was performed or \bold{station} indicating that the time-ordered non-decimated wavelet transform was done.}
\item{date}{The date that the transform was performed or the wd was modified.}
\item{bc}{how the boundaries were handled}
}
\details{
To retain your sanity the C and D coefficients should be extracted by the \code{\link{accessC}} and \code{\link{accessD}} functions and inserted using the \code{\link{putC}} and \code{\link{putD}} functions (or more likely, their methods), rather than by the \code{$} operator.
 
Mind you, if you want to muck about with coefficients directly, then you'll have to do it yourself by working out what the fl.dbase list means (see \code{\link{first.last}} for a description.)
 
Note the \emph{time-ordered non-decimated wavelet transform} used to be called the \emph{stationary wavelet transform}. In fact, the non-decimated transform has several possible names and has been reinvented many times. There are two versions of the non-decimated transform: the coefficients are the same in each version just ordered differently within a resolution level. The two transforms are 

\itemize{
\item{The function \code{\link{wd}}() with an argument \code{type="station"} computes the \emph{time-ordered} non-decimated transform (see Nason and Silverman, 1995) which is useful in time-series applications (see e.g. Nason, von Sachs and Kroisandt, 1998).}
\item{The function \code{\link{wst}}() computes the packets ordered non-decimated transform is useful for curve estimation type applications (see e.g. Coifman and Donoho, 1995). }
}
}
\section{GENERATION}{
This class of objects is returned from the \code{\link{wd}} function to represent a (possibly time-ordered non-decimated) wavelet decomposition of a function. Many other functions return an object of class wd.
}
\section{METHODS}{
The wd class of objects has methods for the following generic functions: \code{\link{plot}}, \code{\link{threshold}}, \code{\link{summary}}, \code{\link{print}}, \code{\link{draw}}.
}
\section{RELAEASE}{Version 3.5.3 Copyright Guy Nason 1994 }
\seealso{
\code{\link{wd}}, \code{\link{wst}}
}
\keyword{classes}
\keyword{smooth}
\author{G P Nason}
