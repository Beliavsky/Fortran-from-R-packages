\name{getdetectionlimit}
\alias{getDetectionlimit}
\title{Gets the detection limit stored in the data set}
\description{
  The detection limit of those values below-detection-limit are 
  stored as negative values in compositional dataset. This function 
  extracts that information.
}
\usage{
getDetectionlimit(x,dl=attr(x,"detectionlimit"))
}
\arguments{
  \item{x}{a data set}
  \item{dl}{a default to replace the information in the dataset}
}
\value{
  a matrix in the same shape as x, with a positive value (the detection limit)
  where available, and NA in the other cells. 
}
\details{
  For a proper treatment of truncated data it would be necessary to know
  the detection limit even for observed data. Unfortunately, there
  is no clear way to encode this information without annoying the user.
}
\references{
  Boogaart, K.G. v.d., R. Tolosana-Delgado, M. Bren (2006) Concepts for
  handling of zeros and missing values in compositional data, in
  E. Pirard (ed.) (2006)Proceedings of the IAMG'2006 Annual Conference
  on "Quantitative Geology from multiple sources", September 2006,
  Liege, Belgium, S07-01, 4pages, ISBN 978-2-9600644-0-7, \url{http://www.stat.boogaart.de/Publications/iamg06_s07_01.pdf}
}
\seealso{
  \code{\link{compositions.missings}},\code{\link{zeroreplace}}
}
\examples{
x <- c(2,-0.5,4,3,-0.5,5,BDLvalue,MARvalue,MNARvalue)
getDetectionlimit(x)

}
\author{K.Gerald van den Boogaart}
\keyword{NA}