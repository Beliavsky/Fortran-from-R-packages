\name{cox.stuart.test}
\alias{cox.stuart.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cox Stuart Trend Test}
%%  ~~function to do ... ~~

\description{
Performs the Cox Stuart test of randomness.
}
\usage{
cox.stuart.test(x, alternative)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector containing the data}
  \item{alternative}{a character string with the alternative hypothesis. Must be one of "\code{two.sided}" (default), "\code{left.sided}" or "\code{right.sided}". You can specify just the initial letter.}
}

\details{
Missing values are removed.

Data is grouped in pairs with the ith observation of the first half paired with the ith observation of the second half of the time-ordered data. If the length of vector X is odd the middle observation is eliminated. The cox stuart test is then simply a sign test applied to these paired data.

The possible values "\code{two.sided}", "\code{left.sided}" and "\code{right.sided}" define the alternative hypothesis.
By using the alternative "\code{two.sided}" the null hypothesis of randomness is tested against either an upward trend or an downward trend.  By using the alternative "\code{left.sided}" the null hypothesis of randomness is tested against an upward trend. By using the alternative "\code{right.sided}" the null hypothesis of randomness is tested against a downward trend. 
}


\value{
A list with class "htest" containing the components:
\item{statistic}{The number of pairs with a signal "+"}
\item{n}{The number of pairs, after eliminanting ties.}
\item{p.value}{the p-value for the test.}
\item{alternative}{a character string describing the alternative hypothesis.}
\item{method}{a character string indicating the test performed.}
\item{data.name}{a character string giving the name of the data.}
}
\references{
Conover, W.J. (1999). \emph{Practical Nonparametric Statistics}, 3rd edition, John Wiley & Sons (p. 166).

Cox, D. R. and Stuart, A. (1955). Some quick sign test for trend
in location and dispersion, \emph{Biometrika}, \bold{42}, 80-95.
%URL: \url{http://www.jstor.org/stable/2333424}

Sprent, P. and Smeeton, N.C. (2007). \emph{Applied Nonparametric Statistical Methods}, 4th ed., Chapman and Hall/CRC Texts in Statistical Science (p. 108).
}
\author{
Ayana Mateus
}
\examples{
##
## Example 1
## Conover (1999)
## The total annual precipitation recorded each year, for 19 years.
##
precipitation <- c(45.25, 45.83, 41.77, 36.26, 45.37, 52.25, 35.37, 57.16, 35.37, 58.32, 
41.05, 33.72, 45.73, 37.90, 41.72, 36.07, 49.83, 36.24, 39.90)
cox.stuart.test(precipitation)

##
## Example 2
## Sweet potato production, harvested in the United States, between 1868 and 1937.
##
data(sweetpotato)
cox.stuart.test(sweetpotato$production)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ cox stuart test }
\keyword{ randomness test }% __ONLY ONE__ keyword per line
