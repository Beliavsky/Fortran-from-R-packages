\name{ccomp}
\alias{ccomp}
\title{Count compositions}
\description{
  A class providing the means to analyse count compositions understood
  as Poisson or multinomial realisation, where the portions are given by
  an unkown Aitchison compositions.
}
\usage{
    ccomp(X,parts=1:NCOL(oneOrDataset(X)),total=NA,warn.na=FALSE,
            detectionlimit=NULL,BDL=NULL,MAR=NULL,MNAR=NULL,SZ=NULL)
          }
\arguments{
  \item{X}{composition or dataset of compositions}
  \item{parts}{vector containing the indices xor names of the columns to be used}
  \item{total}{the total amount to be used, typically 1 or 100}
  \item{warn.na}{should the user be warned in case of NA,NaN or 0
    coding different types of missing values?}
  \item{detectionlimit}{a number, vector or matrix of positive
    numbers giving the detection limit of all values, all columns or
    each value, respectively}
  \item{BDL}{the code for 'Below Detection Limit' in X}
  \item{SZ}{the code for 'Structural Zero' in X}
  \item{MAR}{the code for 'Missing At Random' in X}
  \item{MNAR}{the code for 'Missing Not At Random' in X}
}
\value{
  a vector of class \code{"ccomp"} representing count composition 
  or a matrix of class \code{"ccomp"} representing
  multiple count compositions each in one row.
}

\details{ A count composition contains an indirect observation of an
  Aitchison composition by a Poisson or multinomial variable.  A count
  composition can only contain integer counts. It is assumed that the
  total sum is a an artefact and does not contain information on the
  actual composition. But it does contain information on the precision
  of the relative observation.
  
}
\section{Missing Policy}{
  The policy of treatment of zeroes, missing values and values 
  below detecion limit is explained in depth in \link{compositions.missing}. 
  }

\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}

\seealso{ 
  \code{\link{barplot.ccomp}} \code{\link{ccompMultinomialGOF.test}}
  \code{\link{ccompPoissonGOF.test}} \code{\link{cdt.ccomp}}
  \code{\link{cdtInv.ccomp}}
  \code{\link{fitSameMeanDifferentVarianceModel}}
  \code{\link{groupparts.ccomp}} \code{\link{idt.ccomp}}
  \code{\link{idtInv.ccomp}} \code{\link{is.ccomp}}
  \code{\link{mean.ccomp}} \code{\link{names<-.ccomp}}
  \code{\link{names.ccomp}} \code{\link{plot.ccomp}}
  \code{\link{PoissonGOF.test}} \code{\link{rmultinom.ccomp}}
  \code{\link{rnorm.ccomp}} \code{\link{rpois.ccomp}}
  \code{\link{split.ccomp}} \code{\link{totals.ccomp}}
}
\examples{
data(SimulatedAmounts)
plot(acomp(sa.lognormals))
}
\keyword{classes}
