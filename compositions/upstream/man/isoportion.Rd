\name{isoPortionLines}
\alias{isoPortionLines}
\alias{isoPortionLines.acomp}
\alias{isoPortionLines.rcomp}
\alias{isoProportionLines}
\alias{isoProportionLines.acomp}
\alias{isoProportionLines.rcomp}
\title{Isoportion- and Isoproportion-lines }
\description{
  Add lines of equal portion and proportion to ternary diagrams, to serve as reference axis.
}
\usage{
isoPortionLines(...)
\method{isoPortionLines}{acomp}(by=0.2,at=seq(0,1,by=by),...,
                      parts=1:3,total=1,labs=TRUE,lines=TRUE,unit="")
\method{isoPortionLines}{rcomp}(by=0.2,at=seq(0,1,by=by),...,
                      parts=1:3,total=1,labs=TRUE,lines=TRUE,unit="")
isoProportionLines(...)
\method{isoProportionLines}{acomp}(by=0.2,at=seq(0,1,by=by),...,
                      parts=1:3,labs=TRUE,lines=TRUE) 
\method{isoProportionLines}{rcomp}(by=0.2,at=seq(0,1,by=by),...,
                      parts=1:3,labs=TRUE,lines=TRUE) 
          }
\arguments{
  \item{\dots}{graphical arguments}
  \item{at}{numeric in [0,1]: which portions/proportions should be marked?}
  \item{by}{numeric in (0,1]: steps between protions/proportions}
  \item{parts}{numeric vector subset of \{1,2,3\}: the variables to be
    marked}
  \item{total}{the total amount to be used in labeling}
  \item{labs}{logical: plot the labels?}
  \item{lines}{logical: plot the lines?}
  \item{unit}{mark of the units e.g. "\%"}

}
\details{
  Isoportion lines give lines of the same portion of one of the parts,
  while isoproportion line gives lines of the same ratio between two
  parts. The isoproportion lines are straight lines in both the 
  Aitchison and the real geometries of the simplex, while the isoportion 
  lines are not straight in an Aitchison sense (only in the real one).
  However, note that both types of lines remain straight in the real sense
  when perturbed (von Eynatten et al., 2002).
}
\note{
Currently IsoLines only works with individual plots. This is mainly due to the
fact that I have no idea, what the user interface of this function
should look like for multipanel plots. This includes philosophical
problems with the meaning of isoportions in case of marginal plots. 
  }
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{plot.acomp}}
}
\references{
  von Eynatten, H., V. Pawlowsky-Glahn, and J.J. Egozcue (2002) Understanding 
    perturbation on the simplex: a simple method to better visualize and interpret
    compositional data in ternary diagrams. \emph{Mathematical Geology} \bold{34}, 
    249-257\cr
}
\examples{
data(SimulatedAmounts)
plot(acomp(sa.lognormals))
isoPortionLines()
plot(acomp(sa.lognormals),center=TRUE)
isoPortionLines()
plot(rcomp(sa.lognormals))
isoPortionLines()
plot(acomp(sa.lognormals))
isoProportionLines()
plot(acomp(sa.lognormals),center=TRUE)
isoProportionLines()
plot(rcomp(sa.lognormals))
isoProportionLines()
}
\keyword{multivariate}
