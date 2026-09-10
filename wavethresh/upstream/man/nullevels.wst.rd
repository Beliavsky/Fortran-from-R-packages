\name{nullevels.wst}
\alias{nullevels.wst}
\title{Sets whole resolution levels of coefficients equal to zero in a wst object.}
\description{
Sets whole resolution levels of coefficients equal to zero in a \code{\link{wd}} object.
}
\usage{
\method{nullevels}{wst}(wst, levelstonull, \dots)
}
\arguments{
\item{wst}{An object of class \code{\link{wst}}.} 
\item{levelstonull}{An integer vector specifying which resolution levels of coefficients of \code{\link{wst}} that you wish to set to zero. }
\item{\dots}{any other arguments}
}
\details{
Setting whole resolution levels of coefficients to zero can be very useful. For examples, one can construct a linear smoothing method by setting all coefficients above a particular resolution (the \emph{primary resolution} equal to zero. Also setting particular levels equal to zero can also be useful for removing noise which is specific to a particular resolution level (as long as important signal is not also contained at that level). 

To remove individual coefficients on a systematic basis you probably want to look at the \code{\link{threshold}} function. 
}
\value{
An object of class \code{\link{wst}} where the coefficients in resolution levels specified by \code{levelstonull} have been set to zero. 
}
\section{RELEASE}{Version 3.8.1 Copyright Guy Nason 1997 }
\seealso{
\code{\link{nullevels}}, \code{\link{wst}}, \code{\link{wst.object}}, \code{\link{threshold}}. 
}
\examples{
#
# Look at the examples for \code{\link{nullevels.wd}}.
# The operation is almost identical except that \code{\link{wst}}
# objects are replaced by \code{\link{wd}} ones. 
}
\keyword{manip}
\author{G P Nason}
