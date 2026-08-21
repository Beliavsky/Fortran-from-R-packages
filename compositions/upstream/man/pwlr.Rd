\name{pwlr}
\alias{pwlr}
\alias{pwlrInv}
\title{Pairwise log ratio transform}
\description{
  Compute the pairwise log ratio transform of a (dataset of)
  composition(s), and its inverse.
}
\usage{
      pwlr( x, as.rmult=FALSE, as.data.frame=!as.rmult, \dots)
      pwlrInv( z, orig=gsi.orig(z))
  }
\arguments{
  \item{x}{a composition, not necessarily closed}
  \item{z}{the pwlr-transform of a composition, thus a [D(D-1)/2]-dimensional
    real vector, or a matrix with such many columns}
  \item{as.rmult}{logical; should the output be produced as an rmult object?}
  \item{as.data.frame}{logical; should be as a data.frame? if both are false, rmult will be taken}
  \item{\dots}{currently unused}
  \item{orig}{the original composition, to check consistency and recover component names}
}
\value{
  \code{pwlr} gives the pairwise log ratio transform; accepts a compositional dataset
  \code{pwlrInv} gives a closed composition with the given wplr-transform; accepts a dataset
}
\details{
  The pwlr-transform maps a composition in the $D$-part Aitchison-simplex
  isometrically to a $D(D-1)/2$ dimensonal euclidian vector, computing each possible 
  logratio (accounting for the fact that $log(A/B)=-log(B/A)$, and therefore only one of 
  them is necessary). The data can then
  be analysed in this transformation by multivariate
  analysis tools not relying on the invertibility of the covariance function. 
  The interpretation of
  the results is relatively simple, since each component captures the behaviour of the 
  simple ratio between two party. However redundance between them is extremely high, 
  and any of  \code{\link{alr}}, \code{\link{clr}} or \code{\link{ilr}} transformations 
  may be preferred in most applications.\cr
  
  The pairwise logratio transform is given by
  \deqn{ pwlr(x)_{ij} := \ln\frac{x_i}{x_j}  }{pwlr(x)_{ij} := \ln(x_i/x_j)}.
   
  The inverse requires some explanation, because of the redundance between pwlr scores. 
  Note that for any three components $A,B,C$ it holds that $log(A/C)=log(A/B)+log(B/C)$.
  So, any vector of $D(D-1)/2$ coefficients will not be necessarily a valid pwlr-transformed 
  composition: if these coefficients do not satisfy that kind of relations, the vector is, 
  strictly speaking, not a pwlr and should not be inverted. Nevertheless, the function gives 
  a least-squares inversion, as proposed by Tolosana-Delgado and von Eynatten (2009).
}
\author{R. Tolosana-Delgado \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{clr}},\code{\link{alr}},\code{\link{apt}},
  \url{https://ima.udg.edu/Activitats/CoDaWork03/}
}
\references{
  Aitchison, J. (1986) \emph{The Statistical Analysis of Compositional
    Data} Monographs on Statistics and Applied Probability. Chapman &
  Hall Ltd., London (UK). 416p.
  
  Tolosana-Delgado, R. and H. von Eynatten (2009); Grain-size control 
  on petrographic composition of sediments: compositional regression 
  and rounded zeroes. \emph{Mathematical Geosciences}: 41(8): 869-886. 
  \doi{10.1007/s11004-009-9216-6}.
}
\examples{
(tmp <- pwlr(c(1,2,3)))
pwlrInv(tmp)
}
\keyword{multivariate}
