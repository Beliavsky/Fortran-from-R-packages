\name{maxwell}
\alias{maxwell}
\title{Maxwell's RE coefficient for binary data}
\description{
Computes Maxwell's RE as an index of the interrater agreement of binary data.
}
\usage{
maxwell(ratings)
}
\arguments{
  \item{ratings}{n*2 matrix or dataframe, n subjects 2 raters.}
}
\details{
Missing data are omitted in a listwise way.
}
\value{
A list with class '"irrlist"' containing the following components:
  \item{$method}{a character string describing the method applied for the computation of interrater reliability.}
  \item{$subjects}{the number of subjects examined.}
  \item{$raters}{the number of raters (=2).}
  \item{$irr.name}{a character string specifying the name of the coefficient.}
  \item{$value}{value of RE.}
}
\references{
Maxwell, A.E. (1977). Coefficients of agreement between observers and their interpretation. British Journal of Psychiatry, 130, 79-83.
}
\author{Matthias Gamer}
\seealso{
  \code{\link[irr]{kappa2}}
}
\examples{
data(anxiety)
# Median-split to generate binary data
r1 <- ifelse(anxiety$rater1<median(anxiety$rater1),0,1)
r2 <- ifelse(anxiety$rater2<median(anxiety$rater2),0,1)
maxwell(cbind(r1,r2))
}
\keyword{univar}
