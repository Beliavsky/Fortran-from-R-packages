\name{relInterIntra}
\alias{relInterIntra}
\title{Inter- and intra-rater reliability}
\usage{
 relInterIntra(x, nrater=1, raterLabels=NULL, rho0inter=0.6,
               rho0intra=0.8, conf.level=.95)
}
\arguments{
 \item{x}{Data frame or matrix of rater by object scores}
 \item{nrater}{Number of raters}
 \item{raterLabels}{Labels for the raters or methods}
 \item{rho0inter}{Null hypothesis value for the inter-rater reliability coefficient}
 \item{rho0intra}{Null hypothesis value for the intra-rater reliability coefficient}
 \item{conf.level}{Confidence level for the one-sided confidence interval reported}
}
\description{
 \samp{relInterIntra} calculates inter- and intra-rater reliability coefficients.
}
\value{nil}
\author{Tore Wentzel-Larsen}
\references{
Eliasziw, M., Young, S.L., Woodbury, M.G., & Fryday-Field, K. (1994). Statistical methodology for the concurrent assessment of interrater and intrarater reliability: Using goniometric measurements as an example. Physical Therapy, 74, 777-788.
}
\examples{
# testing code for the Goniometer data from the article:
 table4<-matrix(c(
  -2,16,5,11,7,-7,18,4,0,0,-3,3,7,-6,1,-13,2,4,-10,8,7,-3,-5,5,0,7,-8,1,-3,
  0,16,6,10,8,-8,19,5,-3,0,-2,-1,9,-7,1,-14,1,4,-9,9,6,-2,-5,5,-1,6,-8,1,-3,
  1,15,6,10,6,-8,19,5,-2,-2,-2,1,9,-6,0,-14,0,3,-10,8,7,-4,-7,5,-1,6,-8,2,-3,
  2,12,4,9,5,-9,17,5,-7,1,-4,-1,4,-8,-2,-12,-1,7,-10,2,8,-5,-6,3,-4,4,-10,1,-5,
  1,14,4,7,6,-10,17,5,-6,2,-3,-2,4,-10,-2,-12,0,6,-11,8,7,-5,-8,4,-3,4,-11,-1,-4,
  1,13,4,8,6,-9,17,5,-5,1,-3,1,2,-9,-3,-12,0,4,-10,8,7,-5,-7,4,-4,4,-10,0,-5
  ),ncol=6)
 relInterIntra(x=table4,nrater=2,raterLabels=c('universal','Lamoreux'))
}
\keyword{misc}
