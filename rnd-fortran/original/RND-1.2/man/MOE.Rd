\name{MOE}
\alias{MOE}
\title{Mother of All Extractions}
\description{  \code{MOE} function extracts the risk neutral density based on all models and summarizes the results.}
\usage{
MOE(market.calls, call.strikes, market.puts, put.strikes, call.weights = 1, 
  put.weights = 1, lambda = 1, s0, r, te, y, file.name = "myfile")
}
\arguments{
  \item{market.calls}{market calls (most expensive to cheapest)}
  \item{call.strikes}{strikes for the calls (smallest to largest)}
  \item{market.puts}{market calls (cheapest to most expensive)}
  \item{put.strikes}{strikes for the puts (smallest to largest)}
  \item{call.weights}{Weights for the calls (must be in the same order of calls)}
  \item{put.weights}{Weights for the puts (must be in the same order of puts)}
  \item{lambda}{Penalty parameter to enforce the martingale condition}
  \item{s0}{Current asset value}
  \item{r}{risk free rate}
  \item{te}{time to expiration}
  \item{y}{dividend yield}
  \item{file.name}{File names where analysis is to be saved.  SEE DETAILS!}
}
\details{

The MOE function in a few key strokes extracts the risk neutral density
via various methods and summarizes the results.

This function should only be used for European options.

NOTE:  Three files will be produced: 
       filename will have the pdf version of the results.
       file.namecalls.csv will have the predicted call values.
       file.nameputs.csv will have the predicted put values.
}
\value{
  \item{bsm.mu}{mean of log(S(T)), when S(T) is lognormal}
  \item{bsm.sigma}{SD of log(S(T)), when S(T) is lognormal}
  \item{gb.a}{extracted power parameter, when S(T) is assumed to be a GB rv}
  \item{gb.b}{extracted scale paramter, when S(T) is assumed to be a GB rv}
  \item{gb.v}{extracted first beta paramter, when S(T) is assumed to be a GB rv}
  \item{gb.w}{extracted second beta parameter, when S(T) is assumed to be a GB rv}
  \item{mln.alpha.1}{extracted proportion of the first lognormal.  Second one is 1 - \code{alpha.1} in mixture of lognormals}
  \item{mln.meanlog.1}{extracted mean of the log of the first lognormal in mixture of lognormals}
  \item{mln.meanlog.2}{extracted mean of the log of the second lognormal in mixture of lognormals}
  \item{mln.sdlog.1}{extracted standard deviation of the log of the first lognormal in mixture of lognormals}
  \item{mln.sdlog.2}{extracted standard deviation of the log of the second lognormal in mixture of lognormals}
  \item{ew.sigma }{volatility when using the Edgeworth expansions}
  \item{ew.skew }{normalized skewness when using the Edgeworth expansions}
  \item{ew.kurt }{normalized kurtosis when using the Edgeworth expansions}
  \item{a0}{extracted constant term in the quadratic polynomial of Shimko method}
  \item{a1}{extracted coefficient term of k in the quadratic polynomial of Shimko method}
  \item{a2}{extracted coefficient term of k squared in the quadratic polynomial of Shimko method}
}
\references{

  E. Jondeau and S. Poon and M. Rockinger (2007):
  \emph{Financial Modeling Under Non-Gaussian Distributions} 
  Springer-Verlag, London

}
\author{Kam Hamidieh}
\examples{

###
### You should see that all methods extract the same density!
###

r     = 0.05
te    = 60/365
s0    = 1000
sigma = 0.25
y     = 0.02

strikes     = seq(from = 500, to = 1500, by = 5)
bsm.prices  = price.bsm.option(r =r, te = te, s0 = s0, 
              k = strikes, sigma = sigma, y = y)

calls   = bsm.prices$call
puts    = bsm.prices$put

###
### See where your results will go...
###

getwd()


###
###  Running this may take 1-2 minutes...
###
### MOE(market.calls = calls, call.strikes = strikes, market.puts = puts, 
###    put.strikes = strikes, call.weights = 1, put.weights = 1, 
###    lambda = 1, s0 = s0, r = r, te = te, y = y, file.name = "myfile")
###
### You may get some warning messages.  This happens because the
###    automatic initial value selection sometimes picks values
###    that produce NaNs in the generalized beta density estimation.
###    These messages are often inconsequential.
###
}

