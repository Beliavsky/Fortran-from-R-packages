\name{factorialZ}
\alias{factorialZ}
\alias{chooseZ}
\title{Factorial and Binomial Coefficient as Big Integer}
\description{
  Efficiently compute the factorial \eqn{n!} or a binomial coefficient
  \eqn{{n\choose k}}{choose(n, k)} as big integer (class \code{\link{bigz}}).
}
\usage{
factorialZ(n)
chooseZ(n, k)
}
\arguments{
  \item{n}{non-negative integer (vector), for
    \code{factorialZ}.  For \code{chooseZ}, may be a \code{bigz} big
    integer, also negative.}
  \item{k}{non-negative integer vector.}
}
\value{
  a vector of big integers, i.e., of class \code{\link{bigz}}.
}
\seealso{
  \code{\link{factorial}} and \code{\link{gamma}} in base \R;
}
\examples{
factorialZ(0:10)# 1 1 2 6 ... 3628800
factorialZ(0:40)# larger
factorialZ(200)

n <- 1000
f1000 <- factorialZ(n)
stopifnot(1e-15 > abs(as.numeric(1 - lfactorial(n)/log(f1000))))

system.time(replicate(8, f1e4 <<- factorialZ(10000)))
nchar(as.character(f1e4))# 35660 ... (too many to even look at ..)

chooseZ(1000, 100:102)# vectorizes
chooseZ(as.bigz(2)^120, 10)
n <- c(50,80,100)
k <- c(20,30,40)
## currently with an undesirable warning: % from methods/src/eval.c  _FIXME_
stopifnot(chooseZ(n,k) == factorialZ(n) / (factorialZ(k)*factorialZ(n-k)))
}
\keyword{arith}
