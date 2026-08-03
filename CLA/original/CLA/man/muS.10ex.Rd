\name{muS.10ex}
\alias{muS.10ex}
\docType{data}
\title{10 Assets Example Data from Markowitz & Todd}
\description{
  The simple example Data of Markowitz and Todd (2000); used for
  illustrating the CLA; reused in Bailey and López de Prado (2013).
}
\usage{
data("muS.10ex")
}
\format{
  A list with two components,
  \describe{
    \item{mu}{Named num [1:10] 1.175 1.19 0.396 1.12 0.346 ...\cr
            names : chr [1:10] "X1" "X2" "X3" "X4" ...}
    \item{covar}{num [1:10, 1:10] 0.4076 0.0318 0.0518 0.0566 0.033 ...}
  }
}
\source{
  From \file{http://www.quantresearch.info/CLA_Data.csv.txt} (URL no longer
  working, Aug.2020!) by López de Prado.
}
\references{
  Markowitz, H. M. (1987, 1st ed.) and
  Markowitz, H. M. and Todd, P. G. (2000)
  \emph{Mean-Variance Analysis in Portfolio Choice and Capital Markets},
  page 335.

  Bailey, D. H. and López de Prado, M. (2013)
  An open-source implementation of the critical-line algorithm for portfolio
  optimization, \emph{Algorithms} \bold{6}(1), 169--196;
  \doi{10.3390/a6010169}, p. 16f.
}
\examples{
data(muS.10ex)
str(muS.10ex)

CLA.10ex <- with(muS.10ex, CLA(mu, covar, lB=0, uB=1))
if(require("Matrix"))
  drop0(zapsmall(CLA.10ex$weights_set))
## The results, summarized, as in Bayley and López de Prado (Table 2, p.18) :
with(CLA.10ex, round(cbind(MS_weights[,2:1], lambda=lambdas, t(weights_set)), 3))

CLA.10ex.1c <- with(muS.10ex, CLA(mu, covar, lB=1/100, uB=1))
round(CLA.10ex.1c$weights_set, 3)
}
\keyword{datasets}
