\name{sumprojector}
\alias{sumMissingProjector}
\alias{sumMissingProjector.acomp}
\alias{sumMissingProjector.rcomp}
\alias{sumMissingProjector.aplus}
\alias{sumMissingProjector.rplus}
\alias{sumMissingProjector.rmult}

\title{Compute the global projector to the observed subspace.}
\description{
  Routines to compute the global projector to the observed subspace, 
  down-weighting the subspaces with more missing values.
}
\usage{
sumMissingProjector(x,...)
\method{sumMissingProjector}{acomp}(x,has=is.NMV(x),...)
\method{sumMissingProjector}{aplus}(x,has=is.NMV(x),...)
\method{sumMissingProjector}{rcomp}(x,has=!(is.MAR(x)|is.MNAR(x)),...)
\method{sumMissingProjector}{rplus}(x,has=!(is.MAR(x)|is.MNAR(x)),...)
\method{sumMissingProjector}{rmult}(x,has=is.finite(x),...)
          }
\arguments{
  \item{x}{a dataset of some type containing missings}
  \item{has}{the values to be regarded as non missing}
  \item{\dots}{further generic arguments that might be useful for other
    functions.}
}
\value{
 The matrix of rotation/re-weighting of the original data set, 
 down-weighting the subspaces with more missing values. This matrix is useful
 to obtain estimates of the mean (and variance, in the future) still unbiased 
 in the presence of lost values (only of type MAR, stricly-speaking, but anyway
 useful for any type of missing value, when used with care). This matrix is 
 the Fisher Information in the presence of missing values.
}

\details{
 The function \code{\link{missingProjector}} generates a list of N square 
 matrices of dimension DxD (with N and D respectively 
 equal to the number of rows and columns in \code{x}). Each of these 
 matrices gives the projection of a data row onto its observed sub-space.
 Then, the function \code{sumMissingProjector} takes all these matrices and 
 sums them in a efficient way, generating a "summary" of observed sub-spaces.
}
\section{Missing Policy}{
  No missing policy is given by the routine itself. Its treatment of missing values
depends on the "has" argument.
}

\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}, Raimon
Tolosana-Delgado}

\seealso{
  \code{\link{missingProjector}},
  \code{\link{clr}},\code{\link{rcomp}}, \code{\link{aplus}},
  \code{\link{princomp.acomp}}, 
  \code{\link{plot.acomp}}, \code{\link{boxplot.acomp}},
  \code{\link{barplot.acomp}}, \code{\link{mean.acomp}},
  \code{\link{var.acomp}}, \code{\link{variation.acomp}},
  \code{\link{cov.acomp}}, \code{\link{msd}}
}

\references{
  Boogaart, K.G. v.d., R. Tolosana-Delgado, M. Bren (2006) Concepts for
handling of zeros and missing
  values in compositional data, in E. Pirard (ed.) (2006)Proccedings of
  the IAMG'2006 Annual Conference on "Quantitative Geology from multiple
  sources", September 2006, Liege, Belgium, S07-01, 4pages,
  \url{http://stat.boogaart.de/Publications/iamg06_s07_01.pdf}
  
  
}
\examples{
data(SimulatedAmounts)
sumMissingProjector(acomp(sa.lognormals))
sumMissingProjector(acomp(sa.tnormals))
}
\keyword{multivariate}
