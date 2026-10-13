\name{scores}
\alias{scores}
\title{ Calculate scores of the sample }
\description{

This function calculates normal, t, chi-squared, IQR and MAD scores
of given data.

}
\usage{
scores(x, type = c("z", "t", "chisq", "iqr", "mad"), prob = NA, lim = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{x}{ a vector of data. }
  \item{type}{ "z" calculates normal scores (differences between each value and the mean
divided by sd), "t" calculates t-Student scores (transformed by \code{(z*sqrt(n-2))/sqrt(z-1-t^2)}
formula, "chisq" gives chi-squared scores (squares of differences between values and mean
divided by variance. For the "iqr" type, all values lower than first and greater than third
quartile is considered, and difference between them and nearest quartile divided by IQR are
calculated. For the values between these quartiles, scores are always equal to zero.
"mad" gives differences between each value and median, divided by median absolute deviation.}
  \item{prob}{ If set, the corresponding p-values instead of scores are given. If value is set
to 1, p-value are returned. Otherwise, a logical vector is formed, indicating which values are
exceeding specified probability. In "z" and "mad" types, there is also possibility to set this
value to zero, and then scores are confirmed to (n-1)/sqrt(n) value, according to 
Shiffler (1998). The "iqr" type does not support probabilities, but "lim" value can be specified.}
  \item{lim}{ This value can be set for "iqr" type of scores, to form logical vector, which values
has this limit exceeded. }

}
\value{

A vector of scores, probabilities, or logical vector.

}
\references{ 

Schiffler, R.E (1998). Maximum Z scores and outliers. Am. Stat. 42, 1, 79-80.

}
\author{ Lukasz Komsta, corrections by Alan Richter }
\seealso{ \code{\link{mad}}, \code{\link{IQR}}, \code{\link{grubbs.test}}, }
\examples{
set.seed(1234)
x = rnorm(10)
scores(x)
scores(x,prob=1)
scores(x,prob=0.5)
scores(x,prob=0.1)
scores(x,prob=0.93)
scores(x,type="iqr")
scores(x,type="mad")
scores(x,prob=0)
} 
\keyword{models}
