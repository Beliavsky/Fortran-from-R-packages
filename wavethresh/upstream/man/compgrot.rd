\name{compgrot}
\alias{compgrot}
\title{Compute empirical shift for time ordered non-decimated transforms.}
\description{
Computes the empirical shift required for time-ordered non-decimated transform coefficients to bring them into time order.
}
\usage{
compgrot(J, filter.number, family)
}
\arguments{
\item{J}{The \code{number of levels} in the non-decimated transform where coefficients are to be time-aligned.}
\item{filter.number}{The wavelet filter number to be used, see \code{\link{filter.select}}}
\item{family}{The wavelet family, see \code{\link{filter.select}}}
}
\details{
Time-ordered non-decimated transform coefficients when raw are not in exact time alignment due to the phase of the underlying wavelet. This function returns the shifts that are necessary to apply to each resolution level in the transform to bring back each set of time-ordered coefficients into time alignment. Note that the shifts returned are approximate shifts which work for any Daubechies wavelet. More accurate shifts can be computed using detailed knowledge of the particular wavelet used. 

Each shift is "to the left". I.e. higher indexed coefficients should take the place of lower-indexed coefficients. Periodic boundaries are assumed. 

This realignment is mentioned in Walden and Contreras Cristan, (1997) and Nason, Sapatinas and Sawczenko, (1998). 
}
\value{
A vector containing the shifts that need to be applied to each scale level to return them to the correct time alignment.
 
There are \code{J} entries in the vector. The first entry corresponds to the shift required for the finest level coefficients (i.e. level \code{J-1}) and the last entry corresponds to the coarsest level (i.e. level 0). Entry \code{j} corresponds to the shift required for scale level \code{J-j}. 

}
\note{GROT was the shop started by Reginald Perrin. Unfortunately, GROT stands for "Guy ROTation". }
\section{RELEASE}{Version 3.6 Copyright Guy Nason 1997 }
\seealso{
\code{\link{wst}}, \code{\link{wst.object}}, \code{\link{wpst}}, \code{wpst.object}. 
}
\examples{
#
# Let's see how the resolution levels have to be shifted 
#
compgrot(4, filter.number=10, family="DaubExPhase")
#[1]  2  6 15 31
#
# In other words. Scale level 3 needs to be shifted two units.
# Scale level 2 needs to be shifted 6 units
# Scale level 1 needs to be shifted 15 units
# Scale level 0 needs to be shifted 31 units.
}
\keyword{manip}
\author{G P Nason}
