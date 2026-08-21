\name{scale}
\alias{scale.acomp}
\alias{scale.aplus}
\alias{scale.rcomp}
\alias{scale.rplus}
\alias{scale.rmult}
\alias{scale}
\alias{scale.default}
\title{Normalizing datasets by centering and scaling}
\description{
  The dataset is standardized by optional scaling and centering.
}
\usage{
scale(x, center = TRUE, scale = TRUE,\dots)
\method{scale}{default}(x,center=TRUE, scale=TRUE,\dots)
\method{scale}{acomp}(x,center=TRUE, scale=TRUE,...,robust=getOption("robust"))
\method{scale}{rcomp}(x,center=TRUE, scale=TRUE,...,robust=getOption("robust"))
\method{scale}{aplus}(x,center=TRUE, scale=TRUE,...,robust=getOption("robust"))
\method{scale}{rplus}(x,center=TRUE, scale=TRUE,...,robust=getOption("robust"))
\method{scale}{rmult}(x,center=TRUE, scale=TRUE,...,robust=getOption("robust"))
          }
\arguments{
  \item{x}{a dataset or a single vector of some type}
  \item{center}{logical value or the center to be substracted.}
  \item{scale}{logical value or a scaling factor to for multiplication.}
  \item{robust}{A robustness description. See
    \link{robustnessInCompositions} for details. }
  \item{\dots}{added for generic generality}
}
\value{
  a vector or data matrix, as \code{x} and with the same class, but acordingly transformed.
}
\details{
  scaling is defined in different ways for the different data types. It is
  always performed as an operation in the enclosing vector space. In
  all cases an independent scaling of the different coordinates is not always
  appropriate. This is only done for rplus and rmult geometries. The other three
  geometries are treated with a global scaling, keeping the relative variations
  of every part/amount.

  The scaling factors can be a matrix (for cdt or idt space), a scalar,
  or for the r* geometries vector for scaling the entries
  individually. However scaling the entries individually does not make
  sense in the a* geometries. The operation achieve in the r*-geometries
  is indeed the centering of the a*-geometries. 
  }
  \note{
    Note that the \code{"rcomp"} and \code{"rplus"} objects does not
    preserve their
    geometry during scaling and are therefore reported as \code{"rmult"} objects.
    
    See the documentation in package base for details on
   %\code{\link[base]{scale}} and \code{\link[base]{scale.default}}. 
   \code{base::scale} and \code{base::scale.default}
   These functions are only modified
   to allow the additional robustness parameter. 

  }
  \examples{
  data(SimulatedAmounts)
  plot(scale(acomp(sa.groups)))
  \dontrun{
plot(scale(rcomp(sa.groups)))
}
  plot(scale(aplus(sa.groups)))
  \dontrun{
plot(scale(rplus(sa.groups)))
}
  plot(scale(rmult(sa.groups)))

}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{split}{base}}
}
\keyword{multivariate}
