\name{turning.point.test}
\alias{turning.point.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Turning Point Test
}
\description{
Performs the nonparametric Turning Point test of randomness.  
}
\usage{
turning.point.test(x, alternative)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector containing the data}
  \item{alternative}{a character string specifying the alternative hypothesis. Must be one of "\code{two.sided}" (default), "\code{left.sided}" or "\code{right.sided}".}
}  

\details{
Repeated consecutive observations are removed from data.

The possible values "\code{two.sided}", "\code{left.sided}" and "\code{right.sided}" define the alternative hypothesis.
By using the alternative "\code{two.sided}" the null hypothesis of randomness is tested against either a positive or negative serial correlation between neighbouring observations. }
\value{
A list with class "htest" containing the components:
\item{statistic}{the (normalized) value of the statistic test.}
\item{parameter}{the size of the data, after the remotion of consecutive duplicate values.}
\item{p.value}{the p-value for the test.}
\item{alternative}{a character string describing the alternative hypothesis.}
\item{method}{a character string indicating the test performed.}
\item{data.name}{a character string giving the name of the data.}
\item{tp}{the value of the TP statistic (not shown on screen).}
}
\references{
Brockwell, P.J. and Davis, R.A. (2002). \emph{Introduction to Time Series and Forecasting}, 2nd edition, Springer (p. 36).

Mateus, A. and Caeiro, F. (2013). Comparing several tests of randomness based on the difference of observations. In T. Simos, G. Psihoyios and Ch. Tsitouras (eds.), \emph{AIP Conf. Proc.} \bold{1558}, 809--812.

Moore, G.H. and  Wallis, W.A. (1943). Time Series Significance Tests Based on Signs of Differences. \emph{Journal of the
American Statistical Association}, \bold{38}, 153--154.


%% ~put references to the literature/web site here ~
}
\author{
Ayana Mateus and Frederico Caeiro
}
\examples{
##
## Example 1
##
data(sweetpotato)
turning.point.test(sweetpotato$yield)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ turning point }
\keyword{randomness test }% __ONLY ONE__ keyword per line
