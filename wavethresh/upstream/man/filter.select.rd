\name{filter.select}
\alias{filter.select}
\title{Provide wavelet filter coefficients.}
\description{
This function stores the filter coefficients necessary for doing a discrete wavelet transform (and its inverse), including complex-valued compactly supported
wavelets.
}
\usage{
filter.select(filter.number, family="DaubLeAsymm", constant=1)
}
\arguments{
\item{filter.number}{This selects the desired filter, an integer that takes a value dependent upon the family that you select.
For the complex-valued wavelets in the Lina-Mayrand family, the filter number takes the form x.y where x is the number of vanishing moments (3, 4, or 5) and y is the solution number (1 for x = 3 or 4 vanishing moments; 1, 2, 3, or 4 for x = 5 vanishing moments). Note: this argument has a different meaning
for Littlewood-Paley wavelets, see the note below in the Details section.}

\item{family}{
This selects the basic family that the wavelet comes from. The choices are \bold{DaubExPhase} for Daubechies' extremal phase wavelets, \bold{DaubLeAsymm} for Daubechies' ``least-asymmetric'' wavelets,
\bold{Coiflets} for Coiflets,
\bold{Lawton} for Lawton's complex-valued wavelets (equivalent to Lina-Mayrand 3.1 wavelets), \bold{LittlewoodPaley} for a approximation to Littlewood-Paley wavelets, or \bold{LinaMayrand} for the Lina-Mayrand family of complex-valued Daubechies' wavelets.}
\item{constant}{This constant is applied as a multiplier to all the coefficients. It can be a vector, and so you can adapt the filter coefficients to be whatever you want. (This is feature of negative utility, or ``there is less to this than meets the eye'' as my old PhD supervisor would say [GPN]).} 
}
\details{
This function contains at least three types of filter. Two types can be selected with family set to DaubExPhase, these wavelets are the Haar wavelet (selected by filter.number=1 within this family) and Daubechies ``extremal phase'' wavelets selected by filter.numbers ranging from 2 to 10). Setting family to DaubLeAsymm gives you Daubechies least asymmetric wavelets, but here the filter number ranges from 4 to 10. For Daubechies wavelets, filter.number corresponds to the N of that paper, the wavelets become more regular as the filter.number increases, but they are all of compact support. 

With family equal to ``Coiflets'' the function supports
filter numbers ranging from 1 to 5. Coiflets are wavelets where the
scaling function also has vanishing moments.

With family equal to ``LinaMayrand'', the function returns complex-valued Daubechies wavelets. For odd numbers of vanishing moments, there are symmetric complex-valued wavelets i this family, and for five or more vanishing moments there are multiple distinct complex-valued wavelets, distinguished by their (arbitrary) solution number. At present, Lina-Mayrand wavelets 3.1, 4.1, 5.1, 5.2, 5.3, and 5.4 are available in WaveThresh.

Setting family equal to ``Lawton'' chooses complex-valued wavelets. The only wavelet available is the one with ``filter.number'' equal to 3. 

With family equal to ``LittlewoodPaley'' the Littlewood-Paley wavelet
is used. The scaling function is also the same as (or at least proportional
to, depending on your normalization) that of the Shannon scaling
function, so its an approximation to the Shannon wavelet transform.
The ``filter.number'' argument has a special meaning for the Littlewood-Paley
wavelets: it does not represent vanishing moments here. Instead, it
controls the number of filter taps in the quadrature mirror filter:
typically longer values are better, up to the length of the series.
Increasing it higher than the length of the series does not usually have
much effect. Note: extreme caution should be taken with the
Littlewood-Paley wavelet. This implementation is pure time-domain and
as such can only be thought of as an approximation to a complete
Shannon/LP implementation. For example, in actuality the wavelets
are NOT finite impluse response filters like with Daubechies wavelets.
This means that it is possible for an infinite number of Littlewood Paley
wavelet coefficients to be nonzero. However, computers can not store
an infinite number of coefficients and some will be lost. This is most
noticeable with functions with discontinuities and other homogeneities
but it can also happen with some smooth functions. A way to check how "bad"
is can be is to transform your desired function followed immediately by
the inverse transform and compare the original with the resultant sequence.

The function \code{\link{compare.filters}} can be used to compare two filters. 

}
\value{
Alist is returned with four components describing the filter: 

\item{H}{A vector containing the filter coefficients.} 
\item{G}{A vector containing filter coefficients (if Lawton or Lina-Mayrand wavelets are selected, otherwise this is NULL).} 
\item{name}{A character string containing the name of the filter.} 
\item{family}{A character string containing the family of the filter.} 
\item{filter.number}{The filter number used to select the filter from within a family.} 
}
\note{The (Daubechies) filter coefficients should always sum to sqrt(2). This is a useful check on their validity. }
\section{RELEASE}{Version 3.5.3 Copyright Guy Nason 1994, This version
	originally part of the cthresh release which was merged into
	wavethresh in Oct 2012. Original cthresh version due to Stuart
	Barber }
\seealso{
\code{\link{wd}}, \code{\link{wr}}, \code{\link{wr.wd}}, \code{\link{accessC}}, \code{\link{accessD}}, \code{\link{compare.filters}}, \code{\link{imwd}}, \code{\link{imwr}}, \code{\link{threshold}}, \code{\link{draw}}. 
}
\examples{
#This function is usually called by others.
#However, on occasion you may wish to look at the coefficients themselves.

#
# look at the filter coefficients for N=4 (by default Daubechies'
# least-asymmetric wavelets.)
#
filter.select(4)
#$H:
#[1] -0.07576571 -0.02963553  0.49761867  0.80373875  0.29785780
#[6] -0.09921954 -0.01260397  0.03222310
#
#$G:
#NULL
#
#$name:
#[1] "Daub cmpct on least asymm N=4"
#
#$family:
#[1] "DaubLeAsymm"
#
#$filter.number:
#[1] 4
}
\keyword{utilities}
\author{Stuart Barber and G P Nason}
