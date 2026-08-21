\name{Compound}
\Rdversion{2.1}
\alias{Compound}
\alias{fg.cgb2}
\alias{dl.cgb2}
\alias{pl.cgb2}
\alias{dcgb2}
\alias{pcgb2}
\alias{prcgb2}
\title{
Compound Distribution based on the Generalized Beta Distribution of the Second Kind
}
\description{
Mixture distribution based on the compounding property of the GB2, in short "compound GB2". Decomposition of the GB2 distribution with respect to the left and right tail of the distribution. Calculation of the component densities and cumulative distribution functions. 
Calculation of the compound density function and the compound cumulative distribution function.
}
\usage{
fg.cgb2(x, shape1, scale, shape2, shape3, pl0, decomp="r")
dl.cgb2(x, shape1, scale, shape2, shape3, pl0, decomp="r") 
pl.cgb2(y, shape1, scale, shape2, shape3, pl0, decomp="r", tol=1e-05)
dcgb2(x, shape1, scale, shape2, shape3, pl0, pl, decomp="r")
pcgb2(y, shape1, scale, shape2, shape3, pl0, pl, decomp="r")
prcgb2(y1, y2, shape1, scale, shape2, shape3, pl0, pl, decomp="r", tol=1e-08, 
debug=FALSE)
}
\arguments{
  \item{x}{numeric; can be a vector. The value(s) at which the compound density and the component densities are calculated, \code{x} is positive.}
  \item{y}{numeric; can be a vector. The value(s) at which the compound distribution function and the component distribution functions are calculated.}
  \item{y1, y2}{numeric values.}
  \item{shape1, scale ,shape2, shape3}{numeric; positive parameters of the GB2 distribution.}
  \item{pl0}{numeric; a vector of initial proportions defining the number of components and the weight of each component density in the decomposition. Sums to one.}
  \item{pl}{numeric; a vector of fitted proportions. Sums to one. If \code{pl} is equal to \code{pl0}, we obtain the GB2 distribution.}
  \item{decomp}{string; specifying if the decomposition of the GB2 is done with respect to the right tail ("r") or the left tail ("l") of the distribution.
   By default, \code{decomp} = "r" - right tail decomposition.}
  \item{debug}{logical; By default, \code{debug = FALSE}.}
  \item{tol}{numeric; tolerance with default 0, meaning to iterate until additional terms do not change the partial sum.}
}
\details{
The number of components \eqn{L} is given by the length of the vector \code{pl0}. In our examples \eqn{L=3}. Let \eqn{N} denote the length of the vector \code{x}. Function \code{fg.cgb2} calculates the \eqn{L} gamma factors which multiply the GB2 density in order to obtain the component density \eqn{f_\ell}. These component densities are calculated using the function \code{dl.cgb2}. Function \code{pl.cgb2} calculates the corresponding \eqn{L} cumulative component distribution functions. Function \code{dcgb2} calculates the resulting compound density function. Function \code{pcgb2} calculates the compound cumulative distribution function for a vector of values \code{y} and function \code{prcgb2}, given  2 arguments \code{y1} and \code{y2}, calculates the probability \eqn{P(min(y1,y2) < Y < max(y1,y2))}, where the random variable \eqn{Y} follows a compound GB2 distribution.
}
\value{\code{fg.cgb2} returns a matrix of size \eqn{N \times L} of the Gamma factors, \code{dl.cgb2} returns a matrix of size \eqn{N \times L} of component densities, \code{pl.cgb2} returns a matrix containing the \eqn{L} component cdfs, \code{dcgb2} returns a matrix of size \eqn{N \times 1} of the GB2 compound density function, \code{pcgb2} returns a matrix of size \eqn{N \times 1} of the GB2 compound distribution function and \code{prcgb2} returns a probability between 0 and 1.  
}
\references{
 Graf, M., Nedyalkova, D., Muennich, R., Seger, J. and Zins, S. (2011)
 AMELI Deliverable 2.1: Parametric Estimation of Income Distributions and
 Indicators of Poverty and Social Exclusion.
 \emph{Technical report}, AMELI-Project.
}
\author{
Monique Graf and Desislava Nedyalkova
}
\examples{
#\dontrun{

#\library{cubature}

# GB2 parameters
af <- 5
bf <- 20000
pf <- 0.45 
qf <- 0.75

p0 <- rep(1/3,3)
p1 <- c(0.37,0.43,0.2)

# a vector of values
x <- rep(20000*seq(1,2,length.out=9))

#Gamma components
fg.cgb2(20000,af,bf,pf,qf,p0)
fg.cgb2(Inf,af,bf,pf,qf,p0,"l")

#Component densities
dl.cgb2(x,af,bf,pf,qf,p0)
dl.cgb2(20000,af,bf,pf,qf,p0,"l")

#Component cdf
pl.cgb2(25000,af,bf,pf,qf,p0)

#Compound cdf
pcgb2(x,af,bf,pf,qf,p0,p1)
prcgb2(37000,38000,af,bf,pf,qf,p0,p1,"l")
#}
}
\keyword{
distribution
}


