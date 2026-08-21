\name{colorsForOutliers}
\alias{colorsForOutliers1}
\alias{colorsForOutliers2}
\alias{pchForOutliers1}
\title{Create a color/char palette or for groups of outliers}
\description{
  Conveniance Functions to generate meaningfull color palettes for
  factors representing different types of outliers.
}
\usage{
colorsForOutliers1(outfac, family=rainbow,
                          extreme="cyan",outlier="red",ok="gray40",unknown="blue")
colorsForOutliers2(outfac,use=whichBits(gsi.orSum(levels(outfac))),
                        codes=c(2^outer(c(24,16,8),1:7,"-")),ok="yellow")
pchForOutliers1(outfac,ok='.',outlier='\004',extreme='\003',unknown='\004',...,
  other=c('\001','\002','\026','\027','\010','\011','\012','\013','\014','\015',
    '\016',strsplit("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ","")[[1]])
    )  

          }
\arguments{
  \item{outfac}{a factor given by an OutlierClassifier
    (e.g. \code{\link{OutlierClassifier1}}). \code{colorsForOutliers1}
    is used for the types
    "best","type","outlier","grade". \code{colorsForOutliers2} is used
    for type all. }
  \item{family}{a function generating a color palette from  a numer of
    colors requested.}
  \item{extreme}{The color/char for extrem but not definitivly outlying observations.}
  \item{outlier}{The color/char for detected outliers.}
  \item{unknown}{The color/char for observation with unclear
    classification.}
  \item{other}{The character codes for other outlier classes.}
  \item{ok}{The color/char for nonoutlying usual observations.}
  \item{use}{a numerical vector giving the indices of the bits of the
    output to be represented. The sequence of the bits determins how each
    bit is represented.}
  \item{codes}{The color influences to be used for each bit.}
  \item{\dots}{further codings for other factorlevels}
}
\value{
  a character vector of colors or a numeric vector of plot chars.
}
\details{
  This functions are provided for coveniance to quickly generate a
  palette of reasonable colors or plotting chars for groups of outliers
  classified by OutlierClassifier1.
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{OutlierClassifier1}}
}

\examples{
\dontrun{
data(SimulatedAmounts)
data5 <- acomp(sa.outliers5)
olc <- OutlierClassifier1(data5)
plot(data5,col=colorsForOutliers1(olc)[olc])
olc <- OutlierClassifier1(data5,type="all")
plot(data5,col=colorsForOutliers2(olc)[olc])
}
}
\keyword{multivariate}
