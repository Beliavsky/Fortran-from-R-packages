\name{rcompmargin}
\alias{rcompmargin}
\title{Marginal compositions in real geometry}
\description{
  Compute marginal compositions by amalgamating the rest (additively).
}
\usage{
          rcompmargin(X,d=c(1,2),name="+",pos=length(d)+1,what="data")
          }
\arguments{
  \item{X}{composition or dataset of compositions}
  \item{d}{vector containing the indices xor names of the columns to be kept}
  \item{name}{The new name of the amalgamation column}
  \item{pos}{The position where the new amalgamation column should be
    stored. This defaults to the last column.}
  \item{what}{The role of X either \code{"data"} for data (or means) to be
    transformed or \code{"var"} for variances to be transformed.}
}
\value{
  A closed compositions with class \code{"rcomp"} containing the
  selected variables given by \code{d} and the the amalgamation column.  }
\details{ The amalgamation column is simply computed by adding the 
non-selected components after closing the composition. This is
consistent with the \code{\link{rcomp}} approach and is widely used because
of its easy interpretation. However, it often leads to difficult-to-read
ternary diagrams and is inconsistent with the \code{\link{acomp}}
approach.

With the argument \code{what="var"} the function transformes an rcomp
variance to the resulting variance of the resulting composition.
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}, Raimon
olosana-Delgado}
\seealso{ \code{\link{acompmargin}}, \code{\link{rcomp}}}
\section{Missing Policy}{
    MNAR has the highest priority, MAR next and WZERO (BDL,SZ),- values are
    considered as 0 and reported as BDL in the End.  
}

\references{
  References missing
}
\examples{
data(SimulatedAmounts)
plot.rcomp(sa.tnormals5,margin="rcomp")
plot.rcomp(rcompmargin(sa.tnormals5,c("Cd","Zn")))
plot.rcomp(rcompmargin(sa.tnormals5,c(1,2)))

}
\keyword{multivariate}
