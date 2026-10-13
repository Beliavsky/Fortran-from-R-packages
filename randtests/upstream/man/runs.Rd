\name{Runs}
\alias{druns}
\alias{pruns}
\alias{qruns}
\alias{rruns}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Distribution of the Wald Wolfowitz Runs Statistic
}
\description{
Probability function, distribution function, quantile function and random generation for the distribution of the Runs statistic obtained from samples with \eqn{n_1}{n1} and \eqn{n_2}{n2} elements of each type.
}
\usage{
druns(x, n1, n2, log = FALSE)
pruns(q, n1, n2, lower.tail = TRUE, log.p = FALSE)
qruns(p, n1, n2, lower.tail = TRUE, log.p = FALSE)
rruns(n, n1, n2)
}
\arguments{
  \item{x, q}{a numeric vector of quantiles.}
  \item{p}{a numeric vector of probabilities.}
  \item{n}{number of observations to return.}
  \item{n1, n2}{the number of elements of first and second type, respectively.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P[X \eqn{\le} x] otherwise, P[X > x].}  
}
\details{
The Runs distribution has probability function
\deqn{
P(R=r)=
\left\{
\begin{array}{cc}
\frac{2{n_1-1 \choose r/2-1}{n_2-1 \choose r/2-1}}{{n_1+n_2 \choose n_1}}, &  \mbox{if } r \mbox{ is even}\\
\frac{{n_1-1 \choose (r-1)/2}{n_2-1 \choose (r-3)/2}\,+\,{n_1-1 \choose (r-3)/2}{n_2-1 \choose (r-1)/2}}{{n_1+n_2 \choose n_1}},  &  \mbox{if } r \mbox{ is odd}\\
\end{array}
\right.
%\qquad r=2,3,\ldots, n_1+n_2.
}{P(R=r) = 2 choose(n1-1,r/2-1)choose(n2-1,r/2-1)/choose(n1+n2,n1), if r is even

P(R=r) =
}
for \eqn{r=2,3,\ldots, 2\min(n_1+n_2)+c}{r = 2, 3, \ldots, 2 min(n1+n2)+c} with \eqn{c=0} if \eqn{n_1=n_2}{n1 = n2} or \eqn{c=1} if \eqn{n_1 \neq n_2}{n_1 =! n_2}.


If an element of \code{x} is not integer, the result of \code{druns} is zero.

The quantile is defined as the smallest value \eqn{x} such that \eqn{F(x) \ge p}, where \eqn{F} is the distribution function.
}
\value{
\code{druns} gives the probability function, \code{pruns} gives the distribution function and \code{qruns} gives the quantile function.
%, and rruns generates random deviates.

%The length of the result is determined by nn for rwilcox, and is the maximum of the lengths of the numerical parameters for the other functions.

%The numerical parameters other than nn are recycled to the length of the result. Only the first elements of the logical parameters are used.
}
\references{
Swed, F.S. and Eisenhart, C. (1943). 
Tables for Testing Randomness of Grouping in a Sequence of Alternatives, \emph{Ann. Math Statist.} \bold{14}(1), 66-87.
}
\examples{
##
## Example: Distribution Function
## Creates Table I in Swed and Eisenhart (1943), p. 70,
## with n1 = 2 and n1 <= n2 <= 20
##
m <- NULL
for (i in 2:20){
  m <- rbind(m, pruns(2:5,2,i))  
}
rownames(m)=2:20
colnames(m)=2:5
#
#              2         3         4 5
# 2  0.333333333 0.6666667 1.0000000 1
# 3  0.200000000 0.5000000 0.9000000 1
# 4  0.133333333 0.4000000 0.8000000 1
# 5  0.095238095 0.3333333 0.7142857 1
# 6  0.071428571 0.2857143 0.6428571 1
# 7  0.055555556 0.2500000 0.5833333 1
# 8  0.044444444 0.2222222 0.5333333 1
# 9  0.036363636 0.2000000 0.4909091 1
# 10 0.030303030 0.1818182 0.4545455 1
# 11 0.025641026 0.1666667 0.4230769 1
# 12 0.021978022 0.1538462 0.3956044 1
# 13 0.019047619 0.1428571 0.3714286 1
# 14 0.016666667 0.1333333 0.3500000 1
# 15 0.014705882 0.1250000 0.3308824 1
# 16 0.013071895 0.1176471 0.3137255 1
# 17 0.011695906 0.1111111 0.2982456 1
# 18 0.010526316 0.1052632 0.2842105 1
# 19 0.009523810 0.1000000 0.2714286 1
# 20 0.008658009 0.0952381 0.2597403 1
# 
}
\keyword{ distribution }
%\keyword{ randomness test }% __ONLY ONE__ keyword per line
