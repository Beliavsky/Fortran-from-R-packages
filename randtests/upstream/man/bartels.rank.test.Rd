\name{bartels.rank.test}
\alias{bartels.rank.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bartels Rank Test
}
\description{
Performs the Bartels rank test of randomness.
}
\usage{
bartels.rank.test(x, alternative, pvalue="normal")
}
% alternative = c("two.sided", "less", "greater")
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector containing the observations}
  \item{alternative}{a character string with the alternative hypothesis. Must be one of "\code{two.sided}" (default), "\code{left.sided}" or "\code{right.sided}". You can specify just the initial letter.}
  \item{pvalue}{a character string specifying the method used to compute the p-value. Must be one of \code{normal} (default), \code{beta}, \code{exact} or \code{auto}.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Missing values are removed.

This is the rank version of von Neumann's Ratio Test for Randomness (von Neumann, 1941).

The test statistic RVN is 
\deqn{RVN=\frac{\sum_{i=1}^{n-1}(R_i-R_{i+1})^2}{\sum_{i=1}^{n}\left(R_i-(n+1)/2\right)^2}}{RVN=\sum(R_i-R_{i+1})^2 / \sum(R_i-(n+1)/2)^2}
where \eqn{R_i=rank(X_i), i=1,\dots, n}{R_i=rank(X_i), i=1,...,n}. It is known that \eqn{(RVN-2)/\sigma} is asymptotically standard normal, where \eqn{\sigma^2=\frac{4(n-2)(5n^2-2n-9)}{5n(n+1)(n-1)^2}}{\sigma^2=[4(n-2)(5n^2-2n-9)]/[5n(n+1)(n-1)^2]}. 

The possible \code{alternative} are "\code{two.sided}", "\code{left.sided}" and "\code{right.sided}". By using the alternative "two.sided" the null hypothesis of randomness is tested against nonrandomness. By using the alternative "\code{left.sided}" the null hypothesis of randomness is tested against a trend. By using the alternative "\code{right.sided}" the null hypothesis of randomness is tested against a systematic oscillation.

By default (if \code{pvalue} is not specified), a normal approximation is used to compute the p-value. With \code{beta}, the p-value is computed using an approximation given by the Beta distribution. With \code{exact}, the exact p-value is computed. The option \code{exact} requires the computation of the exact distribution of the statistic test under the null hypothesis and should only be used for small sample sizes (\eqn{n \le 10}).
}
\value{
A list with class "htest" containing the components:
\item{statistic}{the value of the normalized statistic test.}
\item{parameter, n}{the size of the data, after the remotion of consecutive duplicate values.}
\item{p.value}{the p-value of the test.}
\item{alternative}{a character string describing the alternative hypothesis.}
\item{method}{a character string indicating the test performed.}
\item{data.name}{a character string giving the name of the data.}
\item{rvn}{the value of the RVN statistic (not shown on screen).}
\item{nm}{the value of the NM statistic, the numerator of RVN (not shown on screen).}
\item{mu}{the mean value of the RVN statistic (not shown on screen).}
\item{var}{the variance of the RVN statistic (not shown on screen).}
}
\author{Frederico Caeiro}
\references{
Bartels, R. (1982). The Rank Version of von Neumann's Ratio Test for Randomness, \emph{Journal of the American Statistical Association}, \bold{77}(377), 40--46.

Gibbons, J.D. and Chakraborti, S. (2003). \emph{Nonparametric Statistical Inference}, 4th ed. (pp. 97--98).\cr 
URL: \url{http://books.google.pt/books?id=dPhtioXwI9cC&lpg=PA97&ots=ZGaQCmuEUq}

von Neumann, J. (1941). Distribution of the Ratio of the Mean Square Successive Difference to the Variance. \emph{The Annals of Mathematical Statistics} \bold{12}(4), 367--395. doi:10.1214/aoms/1177731677. \url{https://projecteuclid.org/journals/annals-of-mathematical-statistics/volume-12/issue-4/Distribution-of-the-Ratio-of-the-Mean-Square-Successive-Difference/10.1214/aoms/1177731677.full}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{dbartelsrank}}, \code{\link{pbartelsrank}}
}
\examples{
##
## Example 5.1 in Gibbons and Chakraborti (2003), p.98.
## Annual data on total number of tourists to the United States for 1970-1982.
##
years <- 1970:1982
tourists <- c(12362, 12739, 13057, 13955, 14123,  15698, 17523, 18610, 19842, 
      20310, 22500, 23080, 21916)
plot(years, tourists, pch=20)
bartels.rank.test(tourists, alternative="left.sided", pvalue="beta")
# output
#
#  Bartels Ratio Test
#
#data:  tourists 
#statistic = -3.6453, n = 13, p-value = 1.21e-08
#alternative hypothesis: trend 


##
## Example in Bartels (1982).
## Changes in stock levels for 1968-1969 to 1977-1978 (in $A million), deflated by the 
## Australian gross domestic product (GDP) price index (base 1966-1967).
x <- c(528, 348, 264, -20, -167, 575, 410, -4, 430, -122)
bartels.rank.test(x, pvalue="beta")
%GPDprice<-c(83.9, 86.3, 87.3, 88.4, 91.6, 94, 96.9, 100, 102.6, 106.2, 110.9, 117.2, 
%          125.5, 136.3, 155.4, 184.5, 211.8, 234.9, 254.1)
%income<-c(420, 361, 338, 463, 558, 613, 513, 530, 601, 747, 801, 790, 766, 817, 603, 
%          155, 243, 829, 1255)
%x <- income/GPDprice * 100
%bartels.rank.test(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ bartels rank test }
\keyword{ randomness test }% __ONLY ONE__ keyword per line
