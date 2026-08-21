\name{missingsummary}
\alias{missingSummary}
\alias{missingType}
\title{Classify and summarize missing values in a dataset}
\description{
  Routines classifies codes of missing valuesas numbers in objects of the
  compositions package.
}
\usage{
   missingSummary(x,..., vlabs = colnames(x), 
                 mc=attr(x,"missingClassifier"), 
                 values=eval(formals(missingType)$values))
   missingType(x,..., mc=attr(x,"missingClassifier"),
                 values=c("NMV", "BDL", "MAR", "MNAR", "SZ", "Err"))
}
\arguments{
  \item{x}{a dataset which might contain missings}
  \item{\dots}{additional arguments for mc}
  \item{mc}{optionally in missingSummary, an alternate routine to be used 
            instead of \code{missingType} }
  \item{vlabs}{labels for the variables}
  \item{values}{the names of the different types of missings. \code{"Err"} is a
    value that can not be classified e.g. \code{Inf}.}
}
\value{
  \code{missingType} returns a character vector/matrix with the same dimension and
  dimnames as \code{x} giving the type of every value.\cr
  \code{missingSummary} returns a table giving the number of missings of each
  type for each variable.
}
\details{
  The function mainly counts the various types of missing values.
}
\references{
  Boogaart, K.G., R. Tolosana-Delgado, M. Bren (2006) Concepts for the
  handling of zeros and missings in compositional data, \emph{Proceedings of
  IAMG 2006, Liege}

}
\seealso{
  \link{compositions.missing}
}
\examples{
data(SimulatedAmounts)
x <- acomp(sa.lognormals)
xnew <- simulateMissings(x,dl=0.05,MAR=0.05,MNAR=0.05,SZ=0.05)
xnew
missingSummary(xnew)
}
\author{K. Gerald van den Boogaart}
\keyword{NA}
