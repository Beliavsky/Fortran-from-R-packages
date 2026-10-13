\name{robinson}
\alias{robinson}
\title{Robinson's A}
\description{
Computes Robinson's A as an index of the interrater reliability of quantitative data.
}
\usage{
robinson(ratings)
}
\arguments{
  \item{ratings}{n*m matrix or dataframe, n subjects m raters.}
}
\details{
Missing data are omitted in a listwise way.
}
\value{
A list with class '"irrlist"' containing the following components:
  \item{$method}{a character string describing the method applied for the computation of interrater reliability.}
  \item{$subjects}{the number of subjects examined.}
  \item{$raters}{the number of raters.}
  \item{$irr.name}{a character string specifying the name of the coefficient.}
  \item{$value}{coefficient of interrater reliability.}
}
\references{
Robinson, W.S. (1957). The statistical measurement of agreement. American Sociological Review, 22, 17-25.
}
\author{Matthias Gamer}
\seealso{
  \code{\link[irr]{finn}}, 
  \code{\link[irr]{icc}}, 
  \code{\link[irr]{meancor}}
}
\examples{
data(anxiety)
robinson(anxiety)
}
\keyword{univar}
