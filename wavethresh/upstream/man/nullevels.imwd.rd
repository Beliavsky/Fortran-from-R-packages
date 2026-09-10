\name{nullevels.imwd}
\alias{nullevels.imwd}
\title{Sets whole resolution levels of coefficients equal to zero in a imwd object.}
\description{
Sets whole resolution levels of coefficients equal to zero in a \code{\link{imwd.object}}
}
\usage{
\method{nullevels}{imwd}(imwd, levelstonull, \dots)
}
\arguments{
\item{imwd}{An object of class \code{\link{imwd}}.} 
\item{levelstonull}{An integer vector specifying which resolution levels of coefficients of \code{\link{imwd}} that you wish to set to zero. }
\item{\dots}{any other arguments}
}
\details{
Setting whole resolution levels of coefficients to zero can be very useful. For examples, one can construct a linear smoothing method by setting all coefficients above a particular resolution (the \emph{primary resolution} equal to zero. Also setting particular levels equal to zero can also be useful for removing noise which is specific to a particular resolution level (as long as important signal is not also contained at that level). 

Note that this function removes the horiztonal, diagonal and vertical detail coefficients at the resolution level specified. It does not remove the father wavelet coefficients at those resolution levels. 

To remove individual coefficients on a systematic basis you probably want to look at the \code{\link{threshold}} function. 
}
\value{
An object of class \code{\link{imwd}} where the coefficients in resolution levels specified by levelstonull have been set to zero. 
}
\section{RELEASE}{Version 3.9.5 Copyright Guy Nason 1998 }
\seealso{
\code{\link{nullevels}}, \code{\link{imwd}}, \code{\link{imwd.object}}, \code{\link{threshold}}. 
}
\examples{
#
# Do the wavelet transform of the Lennon image
#
data(lennon)
lenimwd <- imwd(lennon)
#
# Set scales (resolution levels) 2, 4 and 6 equal to zero.
#
lenwdNL <- nullevels(lenimwd, levelstonull=c(2,4,6))
#
# Now let's plot the coefficients using a nice blue-heat colour map
#
# You will see that coefficients at levels 2, 4 and 6 are black (i.e. zero)
# You can see that coefficients at other levels are unaffected and still
# show the Lennon coefficients.
#
\dontrun{plot(lenwdNL)}
}
\keyword{manip}

\author{G P Nason}
