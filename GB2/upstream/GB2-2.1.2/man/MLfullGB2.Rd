\name{MLfullGB2}
\Rdversion{2.1}
\alias{MLfullGB2}
\alias{ml.gb2}
\alias{mlh.gb2}
\title{
Maximum Likelihood Estimation of the GB2 Based on the Full Log-likelihood
}
\description{
Performs maximum pseudo-likelihood estimation through the general-purpose optimisation function \code{optim} from package \code{stats}. Two methods of optimization are considered: BFGS and L-BFGS-B 
(see \code{optim} documentation for more details). Initial values of the parameters to be optimized over (\eqn{a}, \eqn{b}, 
\eqn{p} and \eqn{q}) are given from the Fisk distribution and \eqn{p=q=1}. The function to be maximized by \code{optim} is the negative of the full log-likelihood and the gradient is equal to the 
negative of the scores, respectively for the case of a sample of persons and a sample of households. 
}
\usage{
ml.gb2(z, w=rep(1, length(z)), method=1, hess=FALSE) 
mlh.gb2(z, w=rep(1, length(z)), hs=rep(1, length(z)), method=1, hess = FALSE)
}
\arguments{
  \item{z}{numeric; vector of data values.}
  \item{w}{numeric; vector of weights. Must have the same length as \code{z}. By default \code{w} is a vector of 1.}
  \item{hs}{numeric; vector of household sizes. Must have the same length as \code{z}. By default \code{hs} is a vector of 1.}
  \item{method}{numeric; the method to be used by \code{optim}. By default, \code{method = }1 and the used method is BFGS. If \code{method = }2, then the used method is L-BFGS-B.}
  \item{hess}{logical; By default, \code{hess = FALSE}, the hessian matrix is not calculated.}
}
\details{
Function \code{ml.gb2} performs maximum likelihood estimation through the general-purpose optimization function \code{optim} from package \code{stats},
based on the full log-likelihood calculated in a sample of persons. Function \code{mlh.gb2} performs maximum likelihood estimation through the general-purpose optimization function \code{optim} 
from package \code{stats}, based on the full log-likelihood calculated in a sample of households.
}
\value{ 
\code{ml.gb2} and \code{mlh.gb2} return a list with 1 argument: \code{opt1} for the output of the BFGS fit or \code{opt2} for the output of the L-BFGS fit. Further values are given by the values of \code{optim}.
}
\references{
 Graf, M., Nedyalkova, D., Muennich, R., Seger, J. and Zins, S. (2011)
 AMELI Deliverable 2.1: Parametric Estimation of Income Distributions and
 Indicators of Poverty and Social Exclusion.
 \emph{Technical report}, AMELI-Project.
}
\seealso{\code{\link{optim}} for the general-purpose optimization and \code{\link{fisk}} for the Fisk distribution.
}
\author{
Monique Graf
}
\examples{
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

# Fit using the full log-likelihood
fitf <- ml.gb2(inc, w)

# Fitted GB2 parameters
af <- fitf$par[1]
bf <- fitf$par[2]
pf <- fitf$par[3]
qf <- fitf$par[4]

# Likelihood
flik <- fitf$value

# If we want to compare the indicators

# GB2 indicators
indicf <- round(main.gb2(0.6,af,bf,pf,qf), digits=3)
# Empirical indicators
indice <- round(main.emp(inc,w), digits=3)

# Plots
plotsML.gb2(inc,af,bf,pf,qf,w)
}
}
\keyword{distribution}

