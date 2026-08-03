\name{magic-package}
\alias{magic-package}
\docType{package}
\title{
\packageTitle{magic}
}
\description{
\packageDescription{magic}
}
\details{
The DESCRIPTION file:
\packageDESCRIPTION{magic}
\packageIndices{magic}

}
\author{
\packageAuthor{magic}

Maintainer: \packageMaintainer{magic}
}
\references{
  \itemize{
    \item   R. K. S. Hankin 2005.  \dQuote{Recreational mathematics with R:
      introducing the \pkg{magic} package}.  R news, 5(1)
  }
}
\keyword{ package }
\examples{

magic(6)

magicplot(magic(8))

magichypercube.4n(1)

is.alicehypercube(magichypercube.4n(1,d=5),4,give.answers=TRUE)

}
