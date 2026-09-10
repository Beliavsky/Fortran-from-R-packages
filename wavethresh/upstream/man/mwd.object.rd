\name{mwd.object}
\alias{mwd.object}
\title{Multiple wavelet decomposition object (1D) }
\description{
These are objects of class 

\code{mwd}

They represent a decomposition of a function with respect to a multiple wavelet basis. 
}
\details{
To retain your sanity the C and D coefficients should be extracted by the \code{\link{accessC}} and \code{\link{accessD}} functions and put using the \code{\link{putC}} and \code{\link{putD}} functions, rather than by the \code{$} operator. 
}
\value{
The following components must be included in a legitimate `mwd' object. 
\item{C}{a matrix containing each level's smoothed data, each column corresponding to one coefficient vector. The wavelet transform works by applying both a smoothing filter and a bandpass filter to the previous level's smoothed data. The top level contains data at the highest resolution level. Each of these levels are stored one after the other in this matrix. The matrix `\code{fl.dbase$first.last.c}' determines exactly which columns in the matrix, store each level.}
\item{D}{wavelet coefficient matrix. If you were to write down the discrete wavelet transform of a function then columns of D would be the vector coefficients of the wavelet basis function s. Like the C, they are also formed in a pyramidal manner, but stored in a linear matrix. The storage details are to be found in `\code{fl.dbase$first.last.d}'.}
\item{nlevelsWT}{The number of levels in the pyramidal decomposition that produces the coefficients. The precise number of levels depends on the number of different wavelet functions used and the preprocessing method used, as well as the number of data points used.}
\item{fl.dbase}{The first last database associated with this decomposition. This is a list consisting of 2 integers, and 2 matrices. The matrices detail how the coefficients are stored in the C and D components of the `mwd.object'. See the help on \code{\link{mfirst.last}} for more information.}
\item{filter}{a list containing the details of the filter that did the decomposition. See \code{\link{mfilter.select}}.}
\item{type}{either \code{"wavelet"} indicating that the ordinary multiple wavelet transform was performed or \code{"station"} indicating that the non-decimated multiple wavelet transform was done.}
\item{prefilter}{Type of preprocessing or prefilter used. This will be specigic for the type of multiple wavelet used.}
\item{date}{The date that the transform was performed or the mwd object was last modified.}
\item{bc}{how the boundaries were handled} 
}

\section{GENERATION}{
This class of objects is returned from the \code{\link{mwd}} function to represent a multiple wavelet decomposition of a function. Many other functions return an object of class mwd.
}
\section{METHODS}{
The mwd class of objects has methods for the following generic functions: \code{\link{accessC}}, \code{\link{accessD}}, \code{\link{draw}}, \code{\link{plot}}, \code{\link{print}}, \code{\link{putC}}, \code{\link{putD}}, \code{\link{summary}}, \code{\link{threshold}}, \code{\link{wr.mwd}}. 
}
\section{RELEASE}{
Version 3.9.6 (Although Copyright Tim Downie, 1995-6). 
}
\seealso{
\code{\link{accessC.mwd}}, \code{\link{accessD.mwd}}, \code{\link{draw.mwd}}, \code{\link{mfirst.last}}, \code{\link{mfilter.select}}, \code{\link{mwd.object}}, \code{\link{mwr}}, \code{\link{plot.mwd}},\code{\link{print.mwd}}, \code{\link{putC.mwd}}, \code{\link{putD.mwd}}, \code{\link{summary.mwd}}, \code{\link{threshold.mwd}}, \code{\link{wd}}, \code{\link{wr.mwd}}. 
}
\keyword{classes}
\author{Tim Downie}
