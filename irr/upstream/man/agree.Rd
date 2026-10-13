\name{agree}
\alias{agree}
\title{Simple and extended percentage agreement}
\description{
Computes simple and extended percentage agreement among raters.
}
\usage{
agree(ratings, tolerance=0)
}
\arguments{
  \item{ratings}{n*m matrix or dataframe, n subjects m raters.}
  \item{tolerance}{number of successive rating categories that should be regarded as rater agreement (see details).}
}
\details{
Missing data are omitted in a listwise way.\cr
Using extended percentage agreement (tolerance!=0) is only possible for numerical values. If tolerance equals 1, for example, raters differing by one scale degree are interpreted as agreeing.
}
\value{
A list with class '"irrlist"' containing the following components:
  \item{$method}{a character string describing the method applied for the computation of interrater reliability.}
  \item{$subjects}{the number of subjects examined.}
  \item{$raters}{the number of raters.}
  \item{$irr.name}{a character string specifying the name of the coefficient.}
  \item{$value}{coefficient of interrater reliability.}
}
\author{Matthias Gamer}
\seealso{
  \code{\link[irr]{kappa2}}, 
  \code{\link[irr]{kappam.fleiss}}, 
  \code{\link[irr]{kappam.light}} 
}
\examples{
data(video)
agree(video)     # Simple percentage agreement
agree(video, 1)  # Extended percentage agreement
}
\keyword{univar}
