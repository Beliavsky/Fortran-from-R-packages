\name{rcomp}
\alias{rcomp}
\title{Compositions as elements of the simplex embedded in the D-dimensional real space}
\description{
  A class providing a way to analyse compositions in the
  philosophical framework of the Simplex as subset of the \eqn{R^D}.
}
\usage{
          rcomp(X,parts=1:NCOL(oneOrDataset(X)),total=1,warn.na=FALSE,
                detectionlimit=NULL,BDL=NULL,MAR=NULL,MNAR=NULL,SZ=NULL)
          }
\arguments{
  \item{X}{composition or dataset of compositions}
  \item{parts}{vector containing the indices xor names of the columns
    to be used}
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
  a vector of class \code{"rcomp"} representing a closed composition
  or a matrix of class \code{"rcomp"} representing
  multiple closed compositions, by rows.  
}
\details{
  Many multivariate datasets essentially describe amounts of D different
  parts in a whole. This has some important implications justifying to
  regard them as a scale on its own, called a "composition".
  The functions around the class \code{"rcomp"} follow the traditional
  (often statistically inconsistent) approach regarding compositions simply
  as a multivariate vector of positive numbers summing up to 1. This space of
  D positive numbers summing to 1 is traditionally called the D-1-dimensional
  simplex.

  The compositional scale was in-depth analysed by Aitchison
  (1986) and he found serious reasons why compositional data should be
  analysed with a different geometry.  The functions around the class
  \code{"\link{acomp}"} follow his
  approach. However the Aitchison approach based on log-ratios is 
  sometimes criticized (e.g. Rehder and Zier, 2002). It cannot deal with 
  absent parts (i.e. zeros). It is sensitive to large measurement errors
  in small amounts. The Aitchison operations cannot represent simple
  mixture of different compositions. The used transformations
  are not uniformly continuous. Straight lines and ellipses in Aitchison
  space look strangely in ternary diagrams. As all uncritical statistical
  analysis, blind application of logratio-based analysis is sometimes
  misleading. Therefore it is sometimes useful to analyse
  compositional data directly as a multivariate dataset of portions
  summing to 1. However a clear warning must be given that the
  utilisation of almost any kind of
  classical multivariate analysis introduce some kinds of artifacts
  (e.g. Chayes 1960) when applied to compositional data. So, extra care
  and considerable expert knowlegde is needed for the proper
  interpretation of results achieved in this non-Aitchison approach. The
  package tries to lead the user around these artifacts as much as
  possible and gives hints to major pitfalls in the help. However
  meaningless results cannot be fully avoided in this (rather inconsistent) approach. 
  \cr
  A side effect of the procedure is to force the compositions to sum to
  one, which is done by the closure operation \code{\link{clo}} . 
  \cr
  The classes rcomp, acomp, aplus, and rplus are designed in a fashion as similar as
  possible, in order to allow direct comparison between results achieved  
  by the different approaches. Especially the acomp logistic transforms
  \code{\link{clr}}, \code{\link{alr}}, \code{\link{ilr}} are mirrored
  by analogous linear transforms \code{\link{cpt}}, \code{\link{apt}},
  \code{\link{ipt}} in the rcomp class framework. 

}
\section{Missing Policy}{
  Missing and Below Detecion Limit Policy is
  explained in deeper detail in \link{compositions.missing}. 
  }
\author{Raimon Tolosana-Delgado, K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{cpt}}, \code{\link{apt}}, \code{\link{ipt}},
  \code{\link{acomp}}, \code{\link{rplus}},
  \code{\link{princomp.rcomp}},
  \code{\link{plot.rcomp}}, \code{\link{boxplot.rcomp}},
  \code{\link{barplot.rcomp}}, \code{\link{mean.rcomp}},
  \code{\link{var.rcomp}}, \code{\link{variation.rcomp}},
  \code{\link{cov.rcomp}}, \code{\link{msd}},
  \code{\link{convex.rcomp}}, \code{\link{+.rcomp}} 
}

\references{
  Aitchison, J. (1986) \emph{The Statistical Analysis of Compositional
    Data} Monographs on Statistics and Applied Probability. Chapman &
  Hall Ltd., London (UK). 416p.\cr

  Rehder, S. and U. Zier (2001) Letter to the Editor: Comment on 
  ``Logratio Analysis and Compositional Distance'' by J. Aitchison, C. 
  Barcel\'o-Vidal, J.A. Mart\'in-Fern\'a ndez and V. Pawlowsky-Glahn,
  \emph{Mathematical Geology}, \bold{33} (7), 845-848.\cr

  Zier, U. and S. Rehder (2002) Some comments on log-ratio transformation and compositional distance,
  \emph{Terra Nostra}, Schriften der Alfred Wegener-Stiftung, 03/2003\cr
  
    van den Boogaart, K.G. and R. Tolosana-Delgado (2008) "compositions": a unified 
  R package to analyze Compositional Data, \emph{Computers &
    Geosciences}, 34 (4), pages 320-338, \doi{10.1016/j.cageo.2006.11.017}.

  
}
\examples{
data(SimulatedAmounts)
plot(rcomp(sa.tnormals))
}
\keyword{multivariate}
