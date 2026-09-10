\name{nullevels.wd}
\alias{nullevels.wd}
\title{Sets whole resolution levels of coefficients equal to zero in a wd object.}
\description{
Sets whole resolution levels of coefficients equal to zero in a \code{\link{wd.object}}
}
\usage{
\method{nullevels}{wd}(wd, levelstonull, \dots)
}
\arguments{
\item{wd}{An object of class \code{\link{wd}}.} 
\item{levelstonull}{An integer vector specifying which resolution levels of coefficients of \code{\link{wd}} that you wish to set to zero. }
\item{\dots}{any other arguments}
}
\details{
Setting whole resolution levels of coefficients to zero can be very useful. For examples, one can construct a linear smoothing method by setting all coefficients above a particular resolution (the \emph{primary resolution} equal to zero. Also setting particular levels equal to zero can also be useful for removing noise which is specific to a particular resolution level (as long as important signal is not also contained at that level). 

Note that this function removes the horiztonal, diagonal and vertical detail coefficients at the resolution level specified. It does not remove the father wavelet coefficients at those resolution levels. 

To remove individual coefficients on a systematic basis you probably want to look at the \code{\link{threshold}} function. 
}
\value{
An object of class \code{\link{wd}} where the coefficients in resolution levels specified by \code{levelstonull} have been set to zero. 
}
\section{RELEASE}{Version 3.8.1 Copyright Guy Nason 1997 }
\seealso{
\code{\link{nullevels}}, \code{\link{wd}}, \code{\link{wd.object}}, \code{\link{threshold}}. 
}
\examples{
#
# Generate some test data
#
test.data <- example.1()$y
#
# Do wavelet transform of test.data and plot the wavelet coefficients
#
wds <- wd(test.data)
\dontrun{plot(wds)}
#
# Now let us set all the coefficients in ODD resolution levels equal to zero!
#
# This is just to illustrate the capabilities of the function. I cannot
# imagine you wanting to do this in practice!
##
wdsnl <- nullevels(wds, levelstonull = c(1, 3, 5, 7))
#
# Now let's plot the result
#
\dontrun{plot(wdsnl, scaling = "by.level")}
#
# Lo and behold the odd levels have been set to zero!
}
\keyword{manip}
\author{G P Nason}
