\name{gsigetBalStruct}
\alias{gsi.getBalStruct}
\title{An auxiliary functions to compute user-defined ilr and ipt transforms.}
\description{
  Parses the expressions for the balance functions.
}
\usage{
gsi.getBalStruct(descr,names,allowMinus=FALSE,allowOne=FALSE)
}
\arguments{
  \item{descr}{the expression}
  \item{names}{character: the variables in the dataset}
  \item{allowMinus}{whether to allow for \code{-} as a synonym for \code{/}}
  \item{allowOne}{whether to allow \code{1} as a variable name playing the role of a
  a constant.}
}
\value{
  provides a list of pairs giving the components in the Numerator and
  the Denominator of the balances.
}
\details{
  These functions will not be usually called on themselves, but 
  through their wrappers, mainly \code{balanceBase}.
}
\note{It is better not to use gsi.* functions directly since they are internal
  functions of the package. Use their wrappers.}
\author{Raimon Tolosana-Delgado, K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{ilrBase}},\code{\link{ipt}},\code{\link{ilr}},
  \url{https://ima.udg.edu/Activitats/CoDaWork03/}
}
\references{
  Egozcue J.J., V. Pawlowsky-Glahn, G. Mateu-Figueras and
  C. Barcel'o-Vidal (2003) Isometric logratio transformations for
  compositional data analysis. \emph{Mathematical Geology}, \bold{35}(3)
  279-300\cr
  }
\examples{
#gsi.getBalStruct(~A/B/C,c("A","B","C"))

}
\keyword{internal}
