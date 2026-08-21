\name{backtransform}
\alias{backtransform}
\alias{backtransform.rmult}
\alias{gsi.orig}
\alias{gsi.getV}
\title{Automatic common backtransformation for compositions}
\description{
  Functions to automatically determine and compute the relevant back-transformation
  for a rmult object.
}
\usage{
          backtransform(x, as=x)
          backtransform.rmult(x, as=x)
          gsi.orig(x,y=NULL)
          gsi.getV(x,y=NULL)
          }
\arguments{
  \item{x}{an rmult object to be backtransformed; for both \code{gsi.*} functions: an rmult object to extract the relevant information from}
  \item{as}{an rmult object previously obtained with any compositional 
  transformation of this package.}
  \item{y}{for both \code{gsi.*} functions: an alternative object to extract the relevant information
  from, in case that \code{x} does not include it}
}
\value{
  For functions \code{backtransform} or \code{backtransform.rmult}, a corresponding matrix or vector containing the backtransformation of \code{x}. Efforts are taken to keep any extra attributes (beyond, "dim", "dimnames" and "class") the argument "x" may have   \\
  For function \code{gsi.orig}, the original data with a compositional class, if it exists (or NULL otherwise).  \\
  For function \code{gsi.getV}, the transposed, inverse matrix of log-contrasts originally used to forward transform the original composition \code{orig} to its coefficients/coordinates. If it does not exists, the output is NULL.
}
\details{
  The general idea of this package is to analyse the same data with
  different geometric concepts, in a fashion as similar as possible. For each of the
  four concepts there exists a family of transforms expressing the geometry
  in aan appropriate manner. Transformed data can be further analysed, and certain
  results may be back-transformed to the original scale. These functions take 
  care of tracking, constructing and computing the inverse transformation, whichever was the original
  geometry and forward transformation used.

  
}
\author{R. Tolosana-Delgado, K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{cdt}}, \code{\link{idt}}, \code{\link{clr}}, \code{\link{cpt}}, \code{\link{ilt}}, \code{\link{iit}}, \code{\link{ilr}}, \code{\link{ipt}}, \code{\link{alr}}, \code{\link{apt}}
}
\references{
    van den Boogaart, K.G. and R. Tolosana-Delgado (2008) "compositions": a unified 
  R package to analyze Compositional Data, \emph{Computers &
    Geosciences}, 34 (4), pages 320-338, \doi{10.1016/j.cageo.2006.11.017}.
}
\examples{
x <- acomp(1:5)
x
backtransform(ilr(x))
backtransform(clr(x))
backtransform(idt(x))
backtransform(cdt(x))
backtransform(alr(x))
}
\keyword{multivariate}
