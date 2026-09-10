\name{irregwd}
\alias{irregwd}
\title{Irregular wavelet transform (decomposition).}
\description{
This function performs the irregular wavelet transform as described in the paper by Kovac and Silverman. 
}
\usage{irregwd(gd, filter.number=2, family="DaubExPhase", bc="periodic", verbose=FALSE)}
\arguments{
\item{gd}{A grid structure which is the output of the \code{\link{makegrid}} function.}
\item{filter.number}{This selects the smoothness of wavelet that you want to use in the decomposition. By default this is 2, the Daubechies extremal phase orthonormal compactly supported wavelet with 2 vanishing moments.}
\item{family}{specifies the family of wavelets that you want to use. Two popular options are "DaubExPhase" and "DaubLeAsymm" but see the help for \code{\link{filter.select}} for more possibilities.}
\item{bc}{specifies the boundary handling. If \code{bc="periodic"} the default, then the function you decompose is assumed to be periodic on it's interval of definition, if \code{bc="symmetric"} then the function beyond its boundaries is assumed to be a symmetric reflection of the function in the boundary. The symmetric option was the implicit default in releases prior to 2.2.}
\item{verbose}{Controls the printing of "informative" messages whilst the computations progress. Such messages are generally annoying so it is turned off by default.}
}
\details{
If one has irregularly spaced one-dimensional regression data (t,y), say. Then the function \code{\link{makegrid}} interpolates this to a regular grid and then the standard wavelet transform is used to transform the interpolated data. However, unlike the standard wavelet denoising set-up the interpolated data, y, values are correlated. Hence the wavelet coefficients of the interpolated will be correlated (even after using an orthogonal transform). Hence, in particular, the variance of each wavelet coefficient may well be different and so this routine also computes those variances using a fast algorithm (related to the two-dimensional wavelet transform). 

When thresholding with \code{\link{threshold.irregwd}} the threshold function makes use of the information about the variance of each coefficient to modify the variance locally on a coefficient by coefficient basis. 
}
\value{
An object of class \code{\link{irregwd}} which is a list with the following components. 

\item{C}{Vector of sets of successively smoothed versions of the interpolated data (see description of equivalent component of \code{\link{wd.object}} for further information.)}
\item{D}{Vector of sets of wavelet coefficients of the interpolated data at different resolution levels. (see description of equivalent component of \code{\link{wd.object}} for further information.)}
\item{c}{Vector that aids in calculation of variances of wavelet coefficients (used by \code{\link{threshold.irregwd}}).}
\item{nlevelsWT}{The number of resolution levels. This depends on the length of the data vector. If \code{length(data)=2^m}, then there will be m resolution levels. This means there will be m levels of wavelet coefficients (indexed 0,1,2,...,(m-1)), and m+1 levels of smoothed data (indexed 0,1,2,...,m).}
\item{fl.dbase}{There is more information stored in the C and D than is described above. In the decomposition ``extra'' coefficients are generated that help take care of the boundary effects, this database lists where these start and finish, so the "true" data can be extracted.}
\item{filter}{A list containing information about the filter type: Contains the string "wavelet" or "station" depending on which type of transform was performed.}
\item{bc}{How the boundaries were handled.}
\item{date}{The date the transform was performed.}
}

\section{RELESASE}{
3.9.4 Code Copyright Arne Kovac 1997
}
\seealso{
\code{\link{makegrid}}, \code{\link{wd}}, \code{\link{wr.wd}}, \code{\link{accessC}}, \code{\link{accessc}}, \code{\link{accessD}}, \code{\link{putD}}, \code{\link{putC}}, \code{\link{filter.select}}, \code{\link{plot.irregwd}}, \code{\link{threshold.irregwd}}.
}
\examples{
#
# See full examples at the end of the help for makegrid. 
#
}
\keyword{smooth}
\author{Arne Kovac}
