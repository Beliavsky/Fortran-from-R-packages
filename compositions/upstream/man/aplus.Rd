\name{aplus}
\alias{aplus}
\title{Amounts analysed in log-scale}
\description{
  A class to analyse positive amounts in a logistic framework.
}
\usage{
  aplus(X,parts=1:NCOL(oneOrDataset(X)),total=NA,warn.na=FALSE,
          detectionlimit=NULL,BDL=NULL,MAR=NULL,MNAR=NULL,SZ=NULL)
          }
\arguments{
  \item{X}{vector or dataset of positive numbers}
  \item{parts}{vector containing the indices  xor names of the columns to be used}
  \item{total}{a numeric vectors giving the total amounts of each
    dataset. }
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
  a vector of class \code{"aplus"} representing a vector of amounts
  or a matrix of class \code{"aplus"} representing
  multiple vectors of amounts, each vector in one row.  
}
\details{
  Many multivariate datasets essentially describe amounts of D different
  parts in a whole. When the whole is large in relation to the
  considered parts, such that they do not exclude each other, or when
  the total amount of each componenten is indeed determined by the
  phenomenon under investigation and not by sampling artifacts (such as dilution
  or sample preparation), then the parts can be treated as amounts rather
  than as a composition (cf. \code{\link{acomp}}, \code{\link{rcomp}}).
  \cr
  Like compositions, amounts have some important properties. Amounts are
  always positive. An amount of exactly zero essentially means that we have a
  substance of another quality. Different amounts - spanning different
  orders of magnitude  -  are often given in
  different units (ppm, ppb, g/l, vol.\%, mass \%, molar
  fraction). Often, these amounts are also taken as indicators of
  other non-measured components (e.g. K as indicator for potassium feldspar),
  which might be proportional to the measured amount. 
  However, in contrast to compositions, amounts
  themselves do matter. Amounts are typically heavily
  skewed and in many practical cases a log-transform makes their
  distribution roughly symmetric, even normal.
  \cr
  In full analogy to Aitchison's compositions, vector
  space operations are introduced for amounts: the perturbation
  \code{\link{perturbe.aplus}} as a vector space addition (corresponding
  to change of units), the power transformation
  \code{\link{power.aplus}} as scalar multiplication describing the law
  of mass action, and a distance \code{\link{dist}} which is
  independent of the chosen units. The induced vector space is mapped
  isometrically to a classical \eqn{R^D} by a simple log-transformation called
  \code{\link{ilt}}, resembling classical log transform approaches.  
  \cr
  The general approach in analysing aplus objects is thus to perform
  classical multivariate analysis on ilt-transformed coordinates (i.e., logs)
  and to backtransform or display the results in such a way that they
  can be interpreted in terms of the original amounts.   
  \cr
  The class aplus is complemented by the \code{\link{rplus}}, allowing to
  analyse amounts directly as real numbers, and by the classes
  \code{\link{acomp}} and \code{\link{rcomp}} to analyse the same data
  as compositions disregarding the total amounts, focusing on relative
  weights only.
  \cr
  The classes rcomp, acomp, aplus, and rplus are designed as similar as
  possible in order to allow direct comparison between results achieved  
  by the different approaches. Especially the acomp simplex transforms
  \code{\link{clr}}, \code{\link{alr}}, \code{\link{ilr}} are mirrored
  in the aplus class by the single bijective isometric transform \code{\link{ilt}}

}
\seealso{
  \code{\link{ilt}},\code{\link{acomp}}, \code{\link{rplus}},
  \code{\link{princomp.aplus}}, 
  \code{\link{plot.aplus}}, \code{\link{boxplot.aplus}},
  \code{\link{barplot.aplus}}, \code{\link{mean.aplus}},
  \code{\link{var.aplus}}, \code{\link{variation.aplus}},
  \code{\link{cov.aplus}}, \code{\link{msd}}
}
\section{Missing Policy}{
  The policy of treatment of zeroes, missing values and values 
  below detecion limit is explained in depth in \link{compositions.missing}. 
  }
\author{Raimon Tolosana-Delgado, K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\references{
    van den Boogaart, K.G. and R. Tolosana-Delgado (2008) "compositions": a unified 
  R package to analyze Compositional Data, \emph{Computers &
    Geosciences}, 34 (4), pages 320-338, \doi{10.1016/j.cageo.2006.11.017}.
  
}
\examples{
data(SimulatedAmounts)
plot(aplus(sa.lognormals))
}
\keyword{multivariate}
