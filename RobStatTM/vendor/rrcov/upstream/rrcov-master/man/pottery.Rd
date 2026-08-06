\name{pottery}
\alias{pottery}
\alias{pottery.test}
\docType{data}
\title{Archaic Greek Pottery data}
\description{
The Archaic Greek Pottery data set contains data on fragments of Greek
pottery which were classified into two groups according to
their origin: Attic or Eritrean. Six chemical variables, metallic
oxide constituents, were measured: Si, Al, Fe, Ca and Ti. The main data set
consists of 13 Attic objects and 14 Eritrean ones. There is a separate data
set with 13 observations which can be used as a test data set. It consists
of 4 observations classified as "probably Attic" and the remaining 9 as "probably Eritrean".
}
\usage{
data(pottery)
}
\format{
  Two data frames with 27 an 13 observations on the following 7 variables.
  \describe{
    \item{\code{SI}}{Si content}
    \item{\code{AL}}{Al content}
    \item{\code{FE}}{Fe content}
    \item{\code{MG}}{Mg content}
    \item{\code{CA}}{Ca content}
    \item{\code{TI}}{Ti content}
    \item{\code{origin}}{Origin - factor with two levels: \code{Attic} and \code{Eritrean}}
  }
}
\details{
The Archaic Greek Pottery data set was first published by
Stern and Descoeudres (1977) and later reproduced in
Cooper and Weeks (1983) for illustration of linear discriminant analisys.
The data set was used by Pires and Branco (2010) for illustration
of their projection pursuit approach to linear discriminant analysis.

}
\source{
STERN, W. B. and DESCOEUDRES, J.-P. (1977)
X-RAY FLUORESCENCE ANALYSIS OF ARCHAIC GREEK POTTERY
\emph{Archaeometry}, Blackwell Publishing Ltd, \bold{19}, 73--86.
}

\references{
Cooper, R.A. and Weekes, A.J.. 1983
\emph{Data, Models, and Statistical Analysis},
(Lanham, MD: Rowman & Littlefield).

Pires, A. M. and A. Branco, J. (2010)
Projection-pursuit approach to robust linear discriminant analysis
\emph{Journal Multivariate Analysis}, Academic Press, Inc., \bold{101}, 2464--2485.
}

\examples{

data(pottery)
x <- pottery[,c("MG", "CA")]
grp <- pottery$origin

##
## Compute robust location and covariance matrix and
## plot the tolerance ellipses
library(rrcov)
(mcd <- CovMcd(x))
col <- c(3,4)
gcol <- ifelse(grp == "Attic", col[1], col[2])
gpch <- ifelse(grp == "Attic", 16, 1)
plot(mcd, which="tolEllipsePlot", class=TRUE, col=gcol, pch=gpch)

##
## Perform classical LDA and plot the data, 0.975 tolerance ellipses
##  and LDA separation line
##
x <- pottery[,c("MG", "CA")]
grp <- pottery$origin
lda <- LdaClassic(x, grp)
lda
e1 <- getEllipse(loc=lda@center[1,], cov=lda@cov)
e2 <- getEllipse(loc=lda@center[2,], cov=lda@cov)

plot(CA~MG, data=pottery, col=gcol, pch=gpch,
    xlim=c(min(MG,e1[,1], e2[,1]), max(MG,e1[,1], e2[,1])),
    ylim=c(min(CA,e1[,2], e2[,2]), max(CA,e1[,2], e2[,2])))

ab <- lda@ldf[1,] - lda@ldf[2,]
cc <- lda@ldfconst[1] - lda@ldfconst[2]
abline(a=-cc/ab[2], b=-ab[1]/ab[2], col=2, lwd=2)

lines(e1, type="l", col=col[1])
lines(e2, type="l", col=col[2])

##
## Perform robust (MCD) LDA and plot data, classical and
##  robust separation line
##
plot(CA~MG, data=pottery, col=gcol, pch=gpch)
lda <- LdaClassic(x, grp)
ab <- lda@ldf[1,] - lda@ldf[2,]
cc <- lda@ldfconst[1] - lda@ldfconst[2]
abline(a=-cc/ab[2], b=-ab[1]/ab[2], col=2, lwd=2)
abline(a=-cc/ab[2], b=-ab[1]/ab[2], col=4, lwd=2)

rlda <- Linda(x, grp, method="mcd")
rlda
ab <- rlda@ldf[1,] - rlda@ldf[2,]
cc <- rlda@ldfconst[1] - rlda@ldfconst[2]
abline(a=-cc/ab[2], b=-ab[1]/ab[2], col=2, lwd=2)

}
\keyword{datasets}

