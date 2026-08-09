\name{nlfb}
\alias{nlfb}
\title{nlfb: nonlinear least squares modeling by functions}
\usage{
nlfb(
  start,
  resfn,
  jacfn = NULL,
  trace = FALSE,
  lower = -Inf,
  upper = Inf,
  weights = NULL,
  data = NULL,
  ctrlcopy = FALSE,
  control = list(),
  ...
)
}
\arguments{
\item{start}{a numeric vector with all elements present
e.g., start=c(b1=200, b2=50, b3=0.3) 
The start vector for this \code{nlfb}, unlike \code{nlxb}, does
not need to be named.}

\item{resfn}{A function that evaluates the residual vector for 
computing the elements of the sum of squares function at the set of 
parameters \code{start}. Where this function is created by actions on 
a formula or expression in \code{nlxb}, this residual vector will be 
created by evaluation of the 'model - data', rather than 
the conventional 'data - model' approach. The sum of squares is the same.}

\item{jacfn}{A function that evaluates the Jacobian of the sum of squares 
function, that is, the matrix of partial derivatives of the residuals 
with respect to each of the parameters. If NULL (default), uses an 
approximation.

The Jacobian MUST be returned as the attribute "gradient" of this function,
allowing \code{jacfn} to have the same name and be the same code block
as \code{resfn}, which may permit some efficiencies of computation.
}

\item{trace}{TRUE for console output during execution}

\item{lower}{a vector of lower bounds on the parameters. 
If a single number, this will be applied to all. Default \code{-Inf}.}

\item{upper}{a vector of upper bounds on the parameters. If a single number, 
this will be applied to all parameters. Default \code{Inf}.}

\item{weights}{A vector of fixed weights or a function producing one.
See the Details below.}

\item{data}{a data frame of variables used by resfn and jacfn to compute the
required residuals and Jacobian.}

\item{ctrlcopy}{If TRUE use control supplied as is.
           This avoids reprocessing controls.}

\item{control}{a list of control parameters. See nlsr.control().}

\item{...}{additional data needed to evaluate the modeling functions}
}
\value{
{
  A list of the following items: \describe{
  \item{coefficients}{A named vector giving the parameter values at
         the supposed solution.}
  \item{ssquares}{The sum of squared residuals at this set of parameters.}
  \item{resid}{The weighted residual vector at the returned parameters.}
  \item{jacobian}{The jacobian matrix (partial derivatives of residuals w.r.t.
         the parameters) at the returned parameters.}
  \item{feval}{The number of residual evaluations (sum of squares
         computations) used.}
  \item{jeval}{The number of Jacobian evaluations used.}
  \item{weights0}{The weights argument as specified.}
  \item{weights}{The weights vector at the final fit.}
}            
}
}
\description{
A simplified and hopefully robust alternative to finding
the nonlinear least squares minimizer that causes
'formula' to give a minimal residual sum of squares.
}
\details{
nlfb is particularly intended to allow for the
resolution of very ill-conditioned or else near
zero-residual problems for which the regular nls()
function is ill-suited. 

This variant uses a qr solution without forming the sum
of squares and cross products t(J)%*%J

Neither this function nor \code{nlxb} have provision for parameter
scaling (as in the \code{parscale} control of \code{optim} and
package \code{optimx}). This would be more tedious than difficult to
introduce, but does not seem to be a priority feature to add.

The \code{weights} argument can be a vector of fixed weights, in which
case the objective function that will be 
minimized is the sum of squares where each residual is multiplied by the 
square root of the corresponding weight. Default \code{NULL} implies 
unit weights.  \code{weights} may alternatively be a function with header \code{function(parms, resids)} to compute such a vector.
}
\examples{
library(nlsr)
# Scaled Hobbs problem
shobbs.res  <-  function(x){ # scaled Hobbs weeds problem -- residual
  # This variant uses looping
  if(length(x) != 3) stop("shobbs.res -- parameter vector n!=3")
  y  <-  c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
           38.558, 50.156, 62.948, 75.995, 91.972)
  tt  <-  1:12
  res  <-  100.0*x[1]/(1+x[2]*10.*exp(-0.1*x[3]*tt)) - y
}
shobbs.jac  <-  function(x) { # scaled Hobbs weeds problem -- Jacobian
  jj  <-  matrix(0.0, 12, 3)
  tt  <-  1:12
  yy  <-  exp(-0.1*x[3]*tt)
  zz  <-  100.0/(1+10.*x[2]*yy)
  jj[tt,1]   <-   zz
  jj[tt,2]   <-   -0.1*x[1]*zz*zz*yy
  jj[tt,3]   <-   0.01*x[1]*zz*zz*yy*x[2]*tt
  attr(jj, "gradient") <- jj
  jj
}
st <- c(b1=2, b2=1, b3=1) # a default starting vector (named!)
# Default controls, standard Nash-Marquardt algorithm
anlf0 <- nlfb(start=st, resfn=shobbs.res, jacfn=shobbs.jac,
               trace=TRUE, control=list(prtlvl=1))
anlf0

# Hartley with step reduction factor of .2
anlf0h <- nlfb(start=st, resfn=shobbs.res, jacfn=shobbs.jac,
              trace=TRUE, control=list(prtlvl=1, lamda=0, laminc=1.0,
              lamdec=1.0, phi=0, stepredn=0.2))
anlf0h

anlf1bm <- nlfb(start=st, resfn=shobbs.res, jacfn=shobbs.jac, lower=c(2,0,0),
                upper=c(2,6,3), trace=TRUE, control=list(prtlvl=1))
anlf1bm
cat("backtrack using stepredn=0.2\n")
anlf1bmbt <- nlfb(start=st, resfn=shobbs.res, jacfn=shobbs.jac, lower=c(2,0,0),
                upper=c(2,6,3), trace=TRUE, control=list(stepredn=0.2, prtlvl=1))
anlf1bmbt
## Short output
pshort(anlf1bm)
anlf2bm <- nlfb(start=st, resfn=shobbs.res, jacfn=shobbs.jac, lower=c(2,0,0),
                upper=c(2,6,9), trace=TRUE, control=list(prtlvl=1))
anlf2bm
cat("backtrack using stepredn=0.2\n")

anlf2bmbt <- nlfb(start=st, resfn=shobbs.res, jacfn=shobbs.jac, lower=c(2,0,0),
                upper=c(2,6,9), trace=TRUE, control=list(stepredn=0.2, prtlvl=1))
anlf2bmbt
## Short output
pshort(anlf2bm)
}
\author{
J C Nash 2014-7-16   nashjc _at_ uottawa.ca
}
