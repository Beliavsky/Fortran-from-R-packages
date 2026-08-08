\name{pso-package}
\alias{pso-package}
\alias{pso}
\docType{package}
\title{
  A Particle Swarm Optimizer
}
\description{
The package provides an implementation of particle swarm optimization
which is consistent with the standard
PSO 2007 and 2011 by Maurice Clerc et al. Additionally a number of ancillary
routines are provided for easy testing and graphics.
}
\details{
\tabular{ll}{
Package: \tab pso\cr
Type: \tab Package\cr
Version: \tab 1.0.4\cr
Date: \tab 2022-04-12\cr
License: \tab LGPL-3\cr
Depends: \tab methods\cr
}

The core function in the package is \code{\link{psoptim}} which can be
used as a drop in replacement for \code{\link{optim}}. When used without
additional control parameters the implementation is intended to be
equivalent to SPSO 2007 (by M. Clerc et al.).

Control parameters can be specified for SPSO 2011 (in its basic
implementation), to clamp the maximal velocity, provide restarting
when the swarm converges to a region as well as using BFGS as a local
search strategy. See \code{\link{psoptim}} for details.
}
\author{
Maintainer: Claus Bendtsen <papyrus.bendtsen@gmail.com>
}
\keyword{ package }
\keyword{ optimize }
\seealso{
  \code{\link{optim}}.
}
\examples{
\dontrun{
## Some examples of using the functions in the package

## Using basic "optim" interface to minimize a function
set.seed(1)
psoptim(rep(NA,2),function(x) 20+sum(x^2-10*cos(2*pi*x)),
        lower=-5,upper=5,control=list(abstol=1e-8))

## Parabola
p <- test.problem("parabola",10) # one local=global minimum
set.seed(1)
o1 <- psoptim(p,control=list(trace=1,REPORT=50))
show(o1)

set.seed(1)
o2 <- psoptim(p,control=list(trace=1,REPORT=50,w=c(.7,.1)))
show(o2)

set.seed(1)
o3 <- psoptim(p,control=list(trace=1,REPORT=1,hybrid=TRUE))
show(o3) ## hybrid much faster

## Griewank
set.seed(2)
p <- test.problem("griewank",10) # lots of local minima
o1 <- psoptim(p,control=list(trace=1,REPORT=50))
show(o1)

## The above sometimes get stuck in a local minima.
## Adding a restart to increase robustness.
set.seed(2)
o2 <- psoptim(p,control=list(trace=1,REPORT=50,reltol=1e-4))
show(o2)

## An then adding the hybrid
set.seed(2)
o3 <- psoptim(p,control=list(trace=1,REPORT=50,reltol=1e-4,
              hybrid=TRUE,hybrid.control=list(maxit=10)))
show(o3)

## Rosenbrock
set.seed(1)
p <- test.problem("rosenbrock",1)
o1 <- psoptim(p,control=list(trace=1,REPORT=50))
show(o1)

## Change to fully informed
set.seed(1)
o2 <- psoptim(p,control=list(trace=1,REPORT=50,p=1))
show(o2)

## Rastrigin
p <- test.problem("rastrigin",10)
set.seed(1)
o1 <- psoptim(p,control=list(trace=1,REPORT=50))
show(o1)

set.seed(1)
o2 <- psoptim(p,control=list(trace=1,REPORT=50,hybrid=TRUE,
              hybrid.control=list(maxit=10)))
show(o2) # better
plot(o1,xlim=c(0,p@maxf),ylim=c(0,100))
lines(o2,col=2) # and much faster convergence

## Ackley
set.seed(1)
p <- test.problem("ackley",10)
o1 <- psoptim(p,control=list(trace=1,REPORT=50))
show(o1)
}}
