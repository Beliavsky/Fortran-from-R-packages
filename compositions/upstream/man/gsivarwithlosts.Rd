\name{gsi.varwithlosts}
\alias{gsi.varwithlosts}
\title{Internal function: computes variance of compositional data set with missing/zero values}
\description{
  Computes an unbiased estimate of the variance of a compositional data set with some missing and zero values.
}
\usage{
gsi.varwithlosts(x,giveCenter=FALSE)
}
\arguments{
  \item{x}{data set of the classes in \code{\link{compositions}} }
  \item{giveCenter}{a logical. If TRUE the function reports the mean as
    an attribute "center" to the result.}
}
\value{
  A variance matrix using all the observed information, where some of its 
  components have been downweighted to account for the missing values. 
  ATTENTION: function quite slow. It will be called directly by \code{\link{var.acomp}} 
  and similar generic functions only if \code{\link{na.action}} is set to \code{\link{na.pass}}.
}
\note{Do not use gsi.* functions directly since they are internal
  functions of the package.}
\author{R. Tolosana-Delgado, K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\references{
  van den Boogaart, K. G., R. Tolosana-Delgado, and M. Bren (2006) Concepts 
  for handling of zeros and missing values in compositional data, In: 
  \emph{2006 Annual Conference of the International Association for Mathematical 
  Geology (IAMG)}, Universit\'e de Li\`ege, Belgium.

  Bren, M., R. Tolosana-Delgado, and K.G. van den Boogaart (2008) News
  from compositions, the R package. In:  Daunis-i-Estadella, J. and
  Mart\'in-Fern\'a ndez, J.A. (Eds.)
  \emph{Proceedings of the 3rd International Workshop on Compositional Data Analysis},
  Universitat de Girona, Spain. \url{https://ima.udg.edu/Activitats/CoDaWork03/}\cr
}
\examples{
# generates a 3-part compositional data set with missings
A <- matrix(c(0.5,0.3,0.1, 0.3,0.2,0, 0.1, 0, 0.75),nrow=3)
Mvar <- 2*ilrvar2clr(A)
Mcenter <- acomp(c(1,5,2,3))
x <- rnorm.acomp(100,Mcenter,Mvar) 
colnames(x)<-c("A","B","C","D")
# eliminate 20 values completely at random
 slost <- sample(1:length(x), size=30)
 x[slost] <- 0
# compute the variance with losts
#print(vr <- gsi.varwithlosts(cdt(x)))
# compare with Mvar
}
\keyword{internal}
