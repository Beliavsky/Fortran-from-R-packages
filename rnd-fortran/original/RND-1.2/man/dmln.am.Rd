\name{dmln.am}
\alias{dmln.am}
\title{Density of Mixture Lognormal for American Options}
\description{\code{mln.am} is the probability density function of a mixture of three lognormal densities.}
\usage{
dmln.am(x, u.1, u.2, u.3, sigma.1, sigma.2, sigma.3, p.1, p.2)
}
\arguments{
  \item{x}{value at which the denisty is to be evaluated}
  \item{u.1}{log mean of the first lognormal}
  \item{u.2}{log mean of the second lognormal}
  \item{u.3}{log mean of the third lognormal}
  \item{sigma.1}{log standard deviation of the first lognormal}
  \item{sigma.2}{log standard deviation of the second lognormal}
  \item{sigma.3}{log standard deviation of the third lognormal}
  \item{p.1}{weight assigned to the first density}
  \item{p.2}{weight assigned to the second density}
}
\details{
mln is density f(x) = p.1 * f1(x) + p.2 * f2(x) + (1 - p.1 - p.2) * f3(x), where
f1, f2, and f3 are lognormal densities with log means u.1,u.2, and u.3 and
standard deviations sigma.1, sigma.2, and sigma.3 respectively.
}
\value{
  \item{out }{density value at x}}
\references{

  Melick, W. R. and Thomas, C. P. (1997). Recovering an asset's implied pdf from option prices: An
  application to crude oil during the gulf crisis. \emph{Journal of Financial and Quantitative Analysis},
  32(1), 91-115.

}
\author{Kam Hamidieh}
\examples{
###
### Just look at a generic density and see if it integrates to 1.
###

u.1     = 4.2
u.2     = 4.5
u.3     = 4.8
sigma.1 = 0.30
sigma.2 = 0.20
sigma.3 = 0.15
p.1     = 0.25
p.2     = 0.45
x = seq(from = 0, to = 250, by = 0.01)
y = dmln.am(x = x, u.1 = u.1, u.2 = u.2, u.3 = u.3, sigma.1 = sigma.1, sigma.2 = sigma.2, 
            sigma.3 = sigma.3, p.1 = p.1, p.2 = p.2)

plot(y ~ x, type="l")
sum(y * 0.01)

###
### Yes, the sum is near 1.
###
}

