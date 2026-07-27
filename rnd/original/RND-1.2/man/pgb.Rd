\name{pgb}
\alias{pgb}
\title{CDF of Generalized Beta}
\description{
   \code{pgb} is the cumulative distribution function (CDF) of a genaralized beta random variable.
}
\usage{
pgb(x, a, b, v, w)
}
\arguments{
  \item{x}{value at which the CDF is to be evaluated}
  \item{a}{power parameter > 0}
  \item{b}{scale paramter > 0}
  \item{v}{first beta paramter > 0}
  \item{w}{second beta parameter > 0}
}
\details{
  Let B be a beta random variable with parameters v and w.  Then Z = b *(B/(1-B))^(1/a) is a generalized
  beta random variable with parameters (a,b,v,w).
}
\value{
  \item{out }{CDF value at x}
}
\references{
  R.M. Bookstaber and J.B. McDonald (1987)
  A general distribution for describing security price returns.
  \emph{Journal of Business}, 60, 401-424 

  X. Liu and M.B. Shackleton and S.J. Taylor and X. Xu (2007)
  Closed-form transformations from risk-neutral to real-world distributions
  \emph{Journal of Business}, 60, 401-424 

  E. Jondeau and S. Poon and M. Rockinger (2007):
  \emph{Financial Modeling Under Non-Gaussian Distributions} 
  Springer-Verlag, London
}
\author{Kam Hamidieh}
\examples{

#
# What does the cdf of a GB look like?
#

a  = 1
b  = 10
v  = 2
w  = 2

x = seq(from = 0, to = 500, by = 0.01)
y = pgb(x = x, a = a, b = b, v = v, w = w)
plot(y ~ x, type = "l")
abline(h=c(0,1), lty=2) 

}
