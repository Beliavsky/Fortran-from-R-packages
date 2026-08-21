\name{groupparts}
\alias{groupparts}
\alias{groupparts.acomp}
\alias{groupparts.rcomp}
\alias{groupparts.aplus}
\alias{groupparts.rplus}
\alias{groupparts.ccomp}
\title{Group amounts of parts}
\description{
  Groups parts by amalgamation or balancing of their amounts or proportions.
}
\usage{
groupparts(x,\dots)
\method{groupparts}{acomp}(x,\dots,groups=list(...))
\method{groupparts}{rcomp}(x,\dots,groups=list(...))
\method{groupparts}{aplus}(x,\dots,groups=list(...))
\method{groupparts}{rplus}(x,\dots,groups=list(...))
\method{groupparts}{ccomp}(x,\dots,groups=list(...))
          }
\arguments{
  \item{x}{an amount/compositional dataset}
  \item{\dots}{further parameters to use (actually ignored)}
  \item{groups}{a list of numeric xor character vectors, each giving a group of parts} 
  }
\value{
  a new dataset of the same type with each group represented by a single column
}
\details{
  In the real geometry grouping is done by amalgamation (i.e. adding the
  parts). In the Aitchison-geometry grouping is done by taking geometric
  means. The new parts are named by named
  formal arguments. Not-mentioned parts remain ungrouped.
}
\section{Missing Policy}{
  For the real geometries, SZ and BDL are considered as 0, and MAR and MNAR
  are kept as missing of the same type. For the relative geometries, a BDL is a
  special kind of MNAR, whereas a SZ is qualitatively different (thus a
  balance with a SZ has no sense). MAR values transfer their MAR property to the 
  resulting new variable. 
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}, Raimon
Tolosana-Delgado}
\seealso{
  \code{\link{aplus}}
}
\references{
  Egozcue, J.J. and V. Pawlowsky-Glahn (2005) Groups of Parts and their
  Balances in Compositional Data Analysis, Mathematical Geology, in press
}
\examples{
data(SimulatedAmounts)
plot(groupparts(acomp(sa.lognormals5),A=c(1,2),B=c(3,4),C=5))
plot(groupparts(aplus(sa.lognormals5),B=c(3,4),C=5))
plot(groupparts(rcomp(sa.lognormals5),A=c("Cu","Pb"),B=c(2,5)))
hist(groupparts(rplus(sa.lognormals5),1:5))
}
\keyword{multivariate}
