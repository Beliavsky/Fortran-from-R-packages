\name{gsi.recode}
\alias{gsi.recodeM2C}
\alias{gsi.recodeC2M}
\alias{gsi.recodeM2Clean}
\alias{gsi.cleanR}
\title{Internal function: Recode missings with IEEE number and vice versa}
\description{
  The missing codes can be replace with arbitrary value e.g. IEEE
  numbers and viceversa through this interface.
}
\usage{
gsi.recodeM2C(x,y=x,BDL,SZ,MAR,MNAR,NMV)
gsi.recodeC2M(x,y=x,na,nan,ninf,inf,neg,zero,pos)
gsi.recodeM2Clean(x,y=x,BDL=NaN,SZ=NaN,MAR=NaN,MNAR=NA,NMV)
gsi.cleanR(x)
}
\arguments{
  \item{x}{the dataset having missings or IEEE numbers}
  \item{y}{a dataset of similar shape, where the replacment should take
    place}
  \item{BDL}{value to replace for BDL}
  \item{SZ}{value to replace for SZ}
  \item{MAR}{value to replace for MAR}
  \item{MNAR}{value to replace for MNAR}
  \item{NMV}{value to replace for NMV}
  \item{na}{value to replace for \code{NA}}
  \item{nan}{value to replace for \code{NaN}}
  \item{ninf}{value to replace for \code{-Inf}}
  \item{inf}{value to replace for \code{Inf}}
  \item{neg}{value to replace for numbers with \code{x<0}}
  \item{zero}{value to replace for numbers with \code{x==0}}
  \item{pos}{value to replace for numbers with \code{x>0}}
}
\value{
  y with entries replaced. gsi.cleanR replaces all improper numbers
  with 0.
}
\note{Do not use gsi.* functions directly since they are internal
  functions of the package. To use a gsi function, attach the \code{gsi}
environment. }
\details{
  This functions are used internally to transform the different types of
  missings correctly.
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{compositions.missing}}
}
\examples{
#gsi.plain(acomp(c(12,3,4)))
}
\keyword{internal}
