\name{zeroreplace}
\alias{zeroreplace}
\title{Zero-replacement routine}
\description{
  A function to automatically replace rounded zeroes/BDLs in a composition.
}
\usage{
          zeroreplace(x,d=NULL,a=2/3)
          }
\arguments{
  \item{x}{composition or dataset of compositions}
  \item{d}{vector containing the detection limits of each part}
  \item{a}{fraction of the detection limit to be used in replacement}
}
\value{
  an object of the same class as \code{x}, where all WZERO values have been replaced. 
  Output contains a further attribute (named \code{Losts}), 
  with a logical array of the same dimensions as \code{x}, 
  showing which elements were replaced (TRUE) and which were 
  kept unchanged (FALSE).}
\details{
  If \code{d} is given, zeroes from each column of \code{x} are replaced by the 
  corresponding detection limit contained there, scaled
  down by the value of \code{a} (usually a scalar, although if 
  it is a vector it will be recycled with a warning). The variable \code{d} should
  be a vector of length equal to \code{ncol(x)} 
  or a matrix of the same shape as \code{x}.
  \cr
  If \code{d=NULL}, then the detection limit is extracted from the data set,
  if it is available there (i.e., if there are negative numbers). If no negative
  number is present in the data set, and no value is given for \code{d}, the
  result will be equal to \code{x}. See \code{\link{compositions.missings}} for more 
  details on the missing policy.
}
\seealso{
  \code{\link{compositions.missings}},\code{\link{getDetectionlimit}}
}
\references{
  Aitchison, J. (1986) \emph{The Statistical Analysis of Compositional
    Data} Monographs on Statistics and Applied Probability. Chapman &
  Hall Ltd., London (UK). 416p.\cr

  Mart\'in-Fern\'andez, J.A.; Barcel\'o-Vidal, C. and Pawlowsky-Glahn, V. (2003)
  Dealing With Zeros and Missing Values in Compositional Data Sets
  Using Nonparametric Imputation. \emph{Mathematical Geology}, 35 , 253-278

  \url{https://ima.udg.edu/Activitats/CoDaWork03/}\cr

  \url{https://ima.udg.edu/Activitats/CoDaWork05/}

}
\examples{
data(SimulatedAmounts)
x <- acomp(sa.lognormals)
xnew <- simulateMissings(x,dl=0.05,knownlimit=FALSE)
xnew
xrep <- zeroreplace(xnew,0.05)
xrep
}
\keyword{multivariate}
