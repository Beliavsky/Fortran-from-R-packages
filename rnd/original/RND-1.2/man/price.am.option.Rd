\name{price.am.option}
\alias{price.am.option}
\title{Price American Options on Mixtures of Lognormals}
\description{  \code{price.am.option} gives the price of a call and a put option 
 at a set strike when the risk neutral density is a mixture of three lognormals.
}
\usage{
price.am.option(k, r, te, w, u.1, u.2, u.3, sigma.1, sigma.2, sigma.3, p.1, p.2)
}
\arguments{
  \item{k}{Strike}
  \item{r}{risk free rate}
  \item{te}{time to expiration}
  \item{w}{Weight, a number between 0 and 1, to weigh the option price bounds}
  \item{u.1}{log mean of the first lognormal}
  \item{u.2}{log mean of the second lognoral}
  \item{u.3}{log mean of the second lognoral}
  \item{sigma.1}{log sd of the first lognormal}
  \item{sigma.2}{log mean of the second lognormal}
  \item{sigma.3}{log mean of the third lognormal}
  \item{p.1}{weight assigned to the first density}
  \item{p.2}{weight assigned to the second density}
}
\details{
  mln is density f(x) = p.1 * f1(x) + p.2 * f2(x) + (1 - p.1 - p.2) * f3(x), where
  f1, f2, and f3 are lognormal densities with log means u.1,u.2, and u.3 and
  standard deviations sigma.1, sigma.2, and sigma.3 respectively.

  Note: Different weight values, w, need to be assigned to whether the call or put
  is in the money or not.  See equations (7) & (8) of Melick and Thomas paper below.
}
\value{
  \item{call.value}{American call value}
  \item{put.value}{American put value}
  \item{expected.f0}{Expected mean value of asset at expiration}
  \item{prob.f0.gr.k}{Probability asset values is greater than strike}
  \item{prob.f0.ls.k}{Probability asset value is less than strike}
  \item{expected.f0.f0.gr.k}{Expected value of asset given asset exceeds strike}
  \item{expected.f0.f0.ls.k}{Expected value of asset given asset is less than strike}
}
\references{

  Melick, W. R. and Thomas, C. P. (1997). Recovering an asset's implied pdf from option prices: An
  application to crude oil during the gulf crisis. \emph{Journal of Financial and Quantitative Analysis},
  32(1), 91-115.

}
\author{Kam Hamidieh}
\examples{

###
###  Set a few parameters and create some
###  American options.  
###

r       = 0.01
te      = 60/365
w.1     = 0.4
w.2     = 0.25
u.1     = 4.2
u.2     = 4.5
u.3     = 4.8
sigma.1 = 0.30
sigma.2 = 0.20
sigma.3 = 0.15
p.1     = 0.25
p.2     = 0.45
theta   = c(w.1,w.2,u.1,u.2,u.3,sigma.1,sigma.2,sigma.3,p.1,p.2)

p.3 = 1 - p.1 - p.2
p.3
expected.f0   =  sum(c(p.1, p.2, p.3) * exp(c(u.1,u.2,u.3) + 
                    (c(sigma.1, sigma.2, sigma.3)^2)/2) )
expected.f0  
 
strikes = 30:170

market.calls = numeric(length(strikes))
market.puts  = numeric(length(strikes))

for (i in 1:length(strikes))
{

  if ( strikes[i] < expected.f0) {
    market.calls[i] = price.am.option(k = strikes[i], r = r, te = te, w = w.1, u.1 = u.1, 
                      u.2 = u.2, u.3 = u.3, sigma.1 = sigma.1, sigma.2 = sigma.2, 
                      sigma.3 = sigma.3, p.1 = p.1, p.2 = p.2)$call.value

    market.puts[i]  = price.am.option(k = strikes[i], r = r, te = te, w = w.2, u.1 = u.1, 
                      u.2 = u.2, u.3 = u.3, sigma.1 = sigma.1, sigma.2 = sigma.2, 
                      sigma.3 = sigma.3, p.1 = p.1, p.2 = p.2)$put.value 
  }  else {

    market.calls[i] = price.am.option(k = strikes[i], r = r, te = te, w = w.2, u.1 = u.1, 
                      u.2 = u.2, u.3 = u.3, sigma.1 = sigma.1, sigma.2 = sigma.2, 
                      sigma.3 = sigma.3, p.1 = p.1, p.2 = p.2)$call.value

    market.puts[i]  = price.am.option(k = strikes[i], r = r, te = te, w = w.1, u.1 = u.1, 
                      u.2 = u.2, u.3 = u.3, sigma.1 = sigma.1, sigma.2 = sigma.2, 
                      sigma.3 = sigma.3, p.1 = p.1, p.2 = p.2)$put.value 
     }

}

###
### Quickly look at the option values...
###

par(mfrow=c(1,2))
plot(market.calls ~ strikes, type="l")
plot(market.puts  ~ strikes, type="l")
par(mfrow=c(1,1))

}

