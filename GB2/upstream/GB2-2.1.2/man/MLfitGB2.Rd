\name{MLfitGB2}
\Rdversion{2.1}
\alias{MLfitGB2}
\alias{main.emp}
\alias{mlfit.gb2}
\title{
Fitting the GB2 by the Method of Maximum Likelihood Estimation and Comparison of the Fitted Indicators with the Empirical Indicators
}
\description{
The function \code{mlfit.gb2} makes a call to \code{ml.gb2} and \code{profml.gb2}. Estimates of the GB2 parameters are calculated using maximum likelihood estimation based on the full and profile log-likelihoods. Empirical estimates of the set of primary indicators of poverty and social inclusion are calculated using the function \code{main.emp} (see package \code{laeken}) and these estimates
are compared with the indicators calculated with the GB2 fitted parameters using the function \code{\link{main.gb2}}.
}
\usage{
main.emp(z, w=rep(1, length(z)))
mlfit.gb2(z, w=rep(1, length(z))) 
}
\arguments{
  \item{z}{numeric; vector of data values.}
  \item{w}{numeric; vector of weights. Must have the same length as \code{z}. By default \code{w} is a vector of 1.}}
\value{
A list containing three different objects. The first is a data frame with the values of the fitted parameters for the full log-likelihood and the profile log-likelihood, 
the values of the two likelihoods, the values of the GB2 estimates of the indicators and the values of the empirical estimates of the indicators. The second and third objects are, respectively,
the fit using the full log-likelihood and the fit using the profile log-likelihood.
}
\seealso{
\code{\link[stats]{optim}} , \code{\link{ml.gb2}}, \code{\link{profml.gb2}}
}
\author{
Monique Graf and Desislava Nedyalkova
}
\examples{
# An example of using the function mlfit.gb2
# See also the examples of ml.gb2 and mlprof.gb2

\dontrun{
library(laeken)
data(eusilc)

# Income
inc <- as.vector(eusilc$eqIncome)

# Weights
w <- eusilc$rb050

# Data set
d <- data.frame(inc, w)
d <- d[!is.na(d$inc),]

# Truncate at 0
inc <- d$inc[d$inc > 0]
w   <- d$w[d$inc > 0]

# ML fit
m1 <- mlfit.gb2(inc,w)

# GB2 fitted parameters and indicators through maximum likelihood estimation
m1[[1]]
# The fit using the full log-likelihood
m1[[2]]
# The fit using the profile log-likelihood
m1[[3]]

# ML fit, when no weights are avalable 
m2 <- mlfit.gb2(inc)
# Results
m2[[1]]
}
}
\keyword{distribution}

