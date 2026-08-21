\name{missingProjector}
\alias{missingProjector}
\alias{missingProjector.acomp}
\alias{missingProjector.rcomp}
\alias{missingProjector.aplus}
\alias{missingProjector.rplus}
\alias{missingProjector.rmult}
\title{Returns a projector the the observed space in case of missings.}
\description{
  Returns projectors on the observed subspace in the presence of missings.
}
\usage{
missingProjector(x,...,by="s")
\method{missingProjector}{acomp}(x,has=is.NMV(x),...,by="s")
\method{missingProjector}{aplus}(x,has=is.NMV(x),...,by="s")
\method{missingProjector}{rcomp}(x,has=!(is.MAR(x)|is.MNAR(x)),...,by="s")
\method{missingProjector}{rplus}(x,has=!(is.MAR(x)|is.MNAR(x)),...,by="s")
}
\arguments{
  \item{x}{a dataset or object of the given class}
  \item{has}{a boolean matrix of the same size indicating nonmissing
    values}
  \item{\dots}{additional arguments for generic purpose only}
  \item{by}{the name of the dataset dimension on \code{has} for
    tensorial computation with tensorA package}
}
\value{
 A dataset of N square matrices of dimension DxD (with N and D respectively 
 equal to the number of rows and columns in \code{x}). Each of these 
 matrices gives the projection of a data row onto its observed sub-space.
 \cr
 The function \code{\link{sumMissingProjector}} takes all these matrices
and sums
 them, generating a "summary" of observed sub-spaces. This matrix is useful
 to obtain estimates of the mean (and variance, in the future) still unbiased 
 in the presence of lost values (only of type MAR, stricly-speaking, but anyway
 useful for any type of missing value, when used with care).
}
\details{
  See the references for details on that function.
}
\references{
  Boogaart, K.G. v.d. (2006) Concepts for handling of zeros and missing
  values in compositional data, in E. Pirard (ed.) (2006)Proccedings of
  the IAMG'2006 Annual Conference on "Quantitative Geology from multiple
  sources", September 2006, Liege, Belgium, S07-01, 4pages,
  \url{http://stat.boogaart.de/Publications/iamg06_s07_01.pdf}
}
\seealso{
  \code{\link{missingsInCompositions}}
}
\examples{
data(SimulatedAmounts)
x <- acomp(sa.lognormals)
xnew <- simulateMissings(x,dl=0.05,MAR=0.05,MNAR=0.05,SZ=0.05)
xnew
plot(missingSummary(xnew))

missingProjector(acomp(xnew))
missingProjector(rcomp(xnew))
missingProjector(aplus(xnew))
missingProjector(rplus(xnew))

}
\author{K.G.van den Boogaart}
\keyword{NA}