\name{plotmissingsummary}
\alias{plot.missingSummary}
\alias{as.missingSummary}
\title{Plot a Missing Summary}
\description{
  Plots a missing summary as a barplot
}
\usage{
\method{plot}{missingSummary}(x,...,main="Missings",legend.text=TRUE,
                    col=c("gray","lightgray","yellow","red","white","magenta"))
 as.missingSummary(x,...)
}
\arguments{
  \item{x}{a missingSummary table with columns representing different
    types of missing}
  \item{\dots}{further graphical parameters to barplot}
  \item{main}{as in \link[graphics]{barplot}}
  \item{legend.text}{as in barplot}
  \item{col}{as in barplot}
}
\value{
  called for its side effect. The return value is not defined.
}
\details{
  The different types of missings are drawn in quasi-self-understandable 
  colors: normal gray for NMV, and lightgray as for BDL (since they contain
  semi-numeric information), yellow (slight warning) for MAR, red (serious 
  warning) for MNAR, white (because they are non-existing) for SZ, and 
  magenta for the strange case of errors.
}
\references{
    See \link{compositions.missings} for more details.
}
\seealso{
  \code{\link{missingSummary}}
}
\examples{
data(SimulatedAmounts)
x <- acomp(sa.lognormals)
xnew <- simulateMissings(x,dl=0.05,MAR=0.05,MNAR=0.05,SZ=0.05)
xnew
plot(missingSummary(xnew))
}
\author{K.Gerald van den Boogaart}
\keyword{multivariate}
