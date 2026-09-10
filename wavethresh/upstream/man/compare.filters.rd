\name{compare.filters}
\alias{compare.filters}
\title{Compares two filters. }
\description{Compares two filters (such as those returned from \code{\link{filter.select}}). This function returns TRUE is they are the same otherwise returns FALSE. 
}
\usage{compare.filters(f1,f2)
}
\arguments{
\item{f1}{Filter, such as that returned by \code{\link{filter.select}} }
\item{f2}{Filter, such as that returned by \code{\link{filter.select}} }
}
\details{
A very simple function. It only needs to check that the \code{family} and \code{filter.number} components of the filter are the same. 
}
\value{
If \code{f1} and \code{f2} are the same the function returns TRUE, otherwise it returns FALSE. 
}
\section{RELEASE}{Version 3.9 Copyright Guy Nason 1998 }
\seealso{
\code{\link{filter.select}}. 
}
\examples{
#
# Create three filters!
#
filt1 <- filter.select(4, family="DaubExPhase")
filt2 <- filter.select(3, family="DaubExPhase")
filt3 <- filter.select(4, family="DaubLeAsymm")
#
# Now let us see if they are the same...
#
compare.filters(filt1, filt2)
# [1] FALSE
compare.filters(filt1, filt3)
# [1] FALSE
compare.filters(filt2, filt3)
# [1] FALSE
#
# Nope, (what a surprise) they weren't. How about
#
compare.filters(filt1, filt1)
# [1] TRUE
#
# Yes, they were the same!
}
\keyword{manip}
\author{G P Nason}
