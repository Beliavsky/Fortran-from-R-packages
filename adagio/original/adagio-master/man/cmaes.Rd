\name{CMAES}
\alias{pureCMAES}
\title{
  Covariance Matrix Adaptation Evolution Strategy
}
\description{
  The CMA-ES (Covariance Matrix Adaptation Evolution Strategy) is an 
  evolutionary algorithm for difficult non-linear non-convex optimization 
  problems in continuous domain.
  The CMA-ES is typically applied to unconstrained or bounded constraint 
  optimization problems, and search space dimensions between three and fifty.
}
\usage{
pureCMAES(par, fun, lower = NULL, upper = NULL, sigma = 0.5,
                    stopfitness = -Inf, stopeval = 1000*length(par)^2, ...)
}
\arguments{
  \item{par}{objective variables initial point.}
  \item{fun}{objective/target/fitness function.}
  \item{lower,upper}{lower and upper bounds for the parameters.}
  \item{sigma}{coordinate wise standard deviation (step size).}
  \item{stopfitness}{stop if fitness < stopfitness (minimization).}
  \item{stopeval}{stop after stopeval number of function evaluations}
  \item{...}{additional parameters to be passed to the function.}
}
\details{
  The CMA-ES implements a stochastic variable-metric method.
  In the very particular case of a convex-quadratic objective function 
  the covariance matrix adapts to the inverse of the Hessian matrix, up
  to a scalar factor and small random fluctuations. The update equations 
  for mean and covariance matrix maximize a likelihood while resembling 
  an expectation-maximization algorithm.
}
\value{
  Returns a list with components \code{xmin} and \code{fmin}.

  Be patient; for difficult problems or high dimensions the function may run 
  for several minutes; avoid problem dimensions of 30 and more!
}
\references{
  Hansen, N. (2011). The CMA Evolution Strategy: A Tutorial.\cr
  \url{https://arxiv.org/abs/1604.00772}

  Hansen, N., D.V. Arnold, and A. Auger (2013). Evolution Strategies. 
  J. Kacprzyk and W. Pedrycz (Eds.). Handbook of Computational Intelligence,
  Springer-Verlag, 2015.
}
\author{
  Copyright (c) 2003-2010 Nikolas Hansen for Matlab code PURECMAES;
  converted to R by Hans W Borchers.
  (Hansen's homepage: www.cmap.polytechnique.fr/~nikolaus.hansen/)
}
\note{
  There are other implementations of Hansen's CMAES in package `cmaes' 
  (simplified form) and in package `parma' as cmaes() (extended form).
}
\seealso{
  \code{cmaes::cmaes}, \code{parma::cmaes}
}
\examples{
\dontrun{
##  Polynomial minimax approximation of data points
##  (see the Remez algorithm)
n <- 10; m <- 101           # polynomial of degree 10; no. of data points
xi <- seq(-1, 1, length = m)
yi <- 1 / (1 + (5*xi)^2)    # Runge's function

pval <- function(p, x)      # Horner scheme
    outer(x, (length(p) - 1):0, "^") \%*\% p

pfit <- function(x, y, n)   # polynomial fitting of degree n
    qr.solve(outer(x, seq(n, 0), "^"), y)

fn1 <- function(p)           # objective function
    max(abs(pval(p, xi) - yi))

pf <- pfit(xi, yi, 10)      # start with a least-squares fitting
sol1 <- pureCMAES(pf, fn1, rep(-200, 11), rep(200, 11))
zapsmall(sol1$xmin)
# [1]  -50.24826    0.00000  135.85352    0.00000 -134.20107    0.00000
# [7]   59.19315    0.00000  -11.55888    0.00000    0.93453

print(sol1$fmin, digits = 10)
# [1] 0.06546780411

##  Polynomial fitting in the L1 norm
##  (or use LP or IRLS approaches)
fn2 <- function(p)
    sum(abs(pval(p, xi) - yi))

sol2 <- pureCMAES(pf, fn2, rep(-100, 11), rep(100, 11))
zapsmall(sol2$xmin)
# [1] -21.93238   0.00000  62.91083   0.00000 -67.84847   0.00000 
# [7]  34.14398   0.00000  -8.11899   0.00000   0.84533

print(sol2$fmin, digits = 10)
# [1] 3.061810639
}
}

\keyword{ optimize }
