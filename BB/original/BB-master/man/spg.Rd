\name{spg}
\alias{spg}
\title{Large-Scale Optimization}
\description{Spectral projected gradient method for large-scale optimization with simple constraints.}
\usage{
     spg(par, fn, gr=NULL, method=3, lower=-Inf, upper=Inf, 
           project=NULL, projectArgs=NULL, 
	   control=list(), quiet=FALSE, alertConvergence=TRUE, ...)
 }
\arguments{
    \item{par}{A real vector argument to \code{fn}, indicating the initial guess for the optimization of nonlinear 
    objective function \code{fn}.}

    \item{fn}{Nonlinear objective function that is to be optimized. 
      A scalar function that takes a real vector as argument and 
      returns a scalar that is the value of the function at that point 
      (see details).}

    \item{gr}{The gradient of the objective function \code{fn} evaluated at the 
      argument.  This is a vector-function that takes a real 
      vector as argument and returns a real vector of the same length.  
      It defaults to "NULL", which means that gradient is evaluated numerically.  
      Computations are dramatically faster in high-dimensional problems when 
      the exact gradient is provided.  See *Example*.} 

    \item{method}{An integer (1, 2, or 3) specifying which Barzilai-Borwein 
      steplength to use.  The default is 3.  See *Details*.} 

    \item{upper}{An upper bound for box constraints.}

    \item{lower}{An lower bound for box constraints.}

    \item{project}{A projection
      function or character string indicating its name. The projection
      function takes a point in \eqn{R^n}{R^n} and 
      projects it onto a region that defines the constraints of the problem.  
      This is a vector-function that takes a real vector as argument and 
      returns a real vector of the same length. See *Details*. 
      If a projection function is not supplied, arguments \code{lower} and
      \code{upper} will cause the use of an internally defined function that
      enforces the implied box constraints.} 

    \item{projectArgs}{A list with arguments to the \code{project}  
      function. See *Details*. }

    \item{control}{A list of control parameters. See *Details*.}

    \item{quiet}{A logical variable (TRUE/FALSE). If \code{TRUE} warnings 
    and some additional information printing are suppressed.  
    Default is \code{quiet = FALSE}
    Note that \code{quiet} and the \code{control} variable \code{trace}  
    affect different printing, so if \code{trace} is not set to \code{FALSE} 
    there will be considerable printed output.}
 
    \item{alertConvergence}{A logical variable. With the default \code{TRUE}
    a warning is issued if convergence is not obtained. When set to \code{FALSE}
    the warning is suppressed.}

    \item{...}{Additional arguments passed to \code{fn} and \code{gr}.
      (Both must accept any specified arguments, either 
      explicitly or by having a \dots argument, but
      they do not need to use them all.)} 
}
\value{A list with the following components:
\item{par}{Parameters that optimize the nonlinear objective function, 
      if convergence is successful.}

\item{value}{The value of the objective function at termination.}

\item{gradient}{L-infinity norm of the projected gradient of the objective function at termination. If convergence is successful, this should be less than \code{gtol}.} 

\item{fn.reduction}{Reduction in the value of the function from its initial value. This is negative in maximization.}

\item{iter}{Number of iterations taken by the algorithm. The gradient is evaluated once each iteration, so the number of gradient evaluations will also be equal to \code{iter}, plus any evaluations necessary for \code{checkGrad}.}

\item{feval}{Number of times the objective \code{fn} was evaluated.}

\item{convergence}{An integer code indicating type of convergence.  \code{0} 
  indicates successful convergence, in which case the projcted gradient is smaller 
  than \code{pgtol} or the change in objective function is smaller than \code{ftol}.  Error codes are: \code{1} indicates that the maximum limit 
  for iterations \code{maxit} has been reached.  \code{2} indicates that 
  maximum limit on function evals has been exceeded.  \code{3} indicates 
  failure due to error in function evaluation. \code{4} indicates failure due 
  to error in gradient evaluation. \code{5} indicates failure due to error in 
  projection.}
   
\item{message}{A text message explaining which termination criterion was used.}
}

\details{
R adaptation, with significant modifications, by Ravi Varadhan, Johns
Hopkins University (March 25, 2008), from the original FORTRAN code of
Birgin, Martinez, and Raydan (2001). The original is available at the 
TANGO project \url{https://www.ime.usp.br/~egbirgin/tango/downloads.php}

A major modification in our R adaptation of the original FORTRAN code is the 
availability of 3 different options for Barzilai-Borwein (BB) 
steplengths: \code{method = 1} is the BB 
steplength used in Birgin, Martinez and Raydan (2000);  \code{method = 2} is 
the other steplength proposed in Barzilai and Borwein's (1988) original paper.  
Finally, \code{method = 3}, is a new steplength, which was first proposed in 
Varadhan and Roland (2008) for accelerating the EM algorithm.  
In fact, Varadhan and Roland (2008) considered 3 similar steplength schemes in 
their EM acceleration work.  Here, we have chosen \code{method = 3} 
as the "default" method.  This method may be slightly slower than the 
other 2 BB steplength schemes, but it generally exhibited more reliable 
convergence to a better optimum in our experiments.  
We recommend that the user try the other steplength schemes if the default 
method does not perform well in their problem.  

Box constraints can be imposed by vectors \code{lower} and \code{upper}. 
Scalar values for \code{lower} and \code{upper} are expanded to apply to 
all parameters. The default \code{lower} is \code{-Inf} and \code{upper} 
is \code{+Inf}, which imply no constraints.

The \code{project} argument provides a way to implement more general constraints
to be imposed on the parameters in \code{spg}. \code{projectArgs} is passed
to the \code{project} function if one is specified.  The first argument of any \code{project} function should be \code{par} and any other arguments should be passed using its argument \code{projectArgs}.    
To avoid confusion it is suggested that user defined \code{project}
functions should not use arguments  \code{lower} and \code{upper}.

The function \code{\link{projectLinear}} incorporates linear equalities and 
inequalities. This function also provides an example of how other projections
might be implemented.
 
Argument \code{control} is a list specifing any changes to default values of 
algorithm control parameters.  Note that the names of these must be 
specified completely.  Partial matching will not work. 
The list items are as follows:
\describe{
\item{M}{A positive integer, typically between 5-20, that controls the monotonicity of the algorithm.  \code{M=1} would enforce strict monotonicity 
in the reduction of L2-norm of \code{fn}, whereas larger values allow for more non-monotonicity.  Global convergence under non-monotonicity is ensured by 
enforcing the Grippo-Lampariello-Lucidi condition (Grippo et al. 1986) in a non-monotone line-search algorithm.  Values of \code{M} 
between 5 to 20 are generally good.  The default is \code{M = 10}.} 

\item{maxit}{The maximum number of iterations.  The default is \code{maxit = 1500}.}

\item{ftol}{Convergence tolerance on the absolute change in objective function between successive iterations.  
Convergence is declared when the change is less than \code{ftol}.  Default is \code{ftol = 1.e-10}.}

\item{gtol}{Convergence tolerance on the infinity-norm of projected gradient \code{gr} evaluated at the current parameter.  
Convergence is declared when the infinity-norm of projected gradient is less
than \code{gtol}.  Default is \code{gtol = 1.e-05}.}

\item{maxfeval}{Maximum limit on the number of function evaluations.  Default is \code{maxfeval = 10000}.}

\item{maximize}{A logical variable indicating whether the objective function is to be maximized.  Default is \code{maximize = FALSE} indicating
minimization.  For maximization (e.g. log-likelihood maximization in statistical modeling), this may be set to \code{TRUE}.}

\item{trace}{A logical variable (TRUE/FALSE). If \code{TRUE}, information on 
the progress of optimization is printed.  Default is \code{trace = !quiet}.}

\item{triter}{An integer that controls the frequency of tracing 
  when \code{trace=TRUE}. Default is \code{triter=10}, which means that
  the objective \code{fn} and the infinity-norm of its projected gradient are
  printed at every 10-th iteration.}   

\item{eps}{A small positive increment used in the finite-difference 
  approximation of gradient.  Default is 1.e-07.}

\item{checkGrad}{\code{NULL} or a  logical variable \code{TRUE/FALSE} 
indicating whether to 
check the provided analytical gradient against a numerical approximation. 
With the default \code{NULL} the gradient is checked if it is estimated to take
less than about ten seconds. A warning will be issued in the case it takes 
longer. The default can be overridden by specifying \code{TRUE} or \code{FALSE}.
It is recommended that this be set to FALSE for high-dimensional problems, 
after making sure that the gradient is correctly specified, possibly by running
once with \code{TRUE} specified.}

\item{checkGrad.tol}{A small positive value use to compare the maximum relative
  difference between a user supplied gradient gr and the numerical approximation
  calculated by grad from package \pkg{numDeriv}. The default is 1.e-06.
  If this value is exceeded then an error message is issued, as it is a 
  reasonable indication of a problem with the user supplied gr. The user can
  either fix the gr function, remove it so the finite-difference 
  approximation is used, or increase the tolerance so the check passes.}
}}

\references{   
    Birgin EG, Martinez JM, and Raydan M (2000): Nonmonotone spectral projected gradient methods on convex sets, \emph{SIAM J Optimization}, 10, 1196-1211.
    
    Birgin EG, Martinez JM, and Raydan M (2001): SPG: software for convex-constrained optimization, \emph{ACM Transactions on Mathematical Software}.
    
   L Grippo, F Lampariello, and S Lucidi (1986),  A nonmonotone line search technique for Newton's method,  \emph{SIAM J on Numerical Analysis}, 23, 707-716. 

   M Raydan (1997),  Barzilai-Borwein gradient method for large-scale unconstrained minimization problem, \emph{SIAM J of Optimization}, 7, 26-33.

   R Varadhan and C Roland (2008), Simple and globally-convergent methods for accelerating the convergence of any EM algorithm,  \emph{Scandinavian J Statistics}, doi: 10.1111/j.1467-9469.2007.00585.x.  

   R Varadhan and PD Gilbert (2009),  BB: An R Package for Solving a Large System of Nonlinear Equations and for Optimizing a High-Dimensional Nonlinear Objective Function, \emph{J. Statistical Software}, 32:4, \url{https://www.jstatsoft.org/v32/i04/}
}
\seealso{
  \code{\link{projectLinear}},
  \code{\link{BBoptim}},
  \code{\link{optim}},
  \code{\link{nlm}},
  \code{\link{sane}},
  \code{\link{dfsane}},
  \code{\link[numDeriv]{grad}}
   }
   
\examples{
sc2.f <- function(x){ sum((1:length(x)) * (exp(x) - x)) / 10}

sc2.g <- function(x){ (1:length(x)) * (exp(x) - 1) / 10}

p0 <- rnorm(50)
ans.spg1 <- spg(par=p0, fn=sc2.f)  # Default is method=3
ans.spg2 <- spg(par=p0, fn=sc2.f, method=1)
ans.spg3 <- spg(par=p0, fn=sc2.f, method=2)
ans.cg <- optim(par=p0, fn=sc2.f, method="CG")  #Uses conjugate-gradient method in "optim"
ans.lbfgs <- optim(par=p0, fn=sc2.f, method="L-BFGS-B")  #Uses low-memory BFGS method in "optim"

# Now we use exact gradient.  
# Computation is much faster compared to when using numerical gradient.
ans.spg1 <- spg(par=p0, fn=sc2.f, gr=sc2.g)

############
# Another example illustrating use of additional parameters to objective function 
valley.f <- function(x, cons) {
  n <- length(x)
  f <- rep(NA, n)
  j <- 3 * (1:(n/3))
  jm2 <- j - 2
  jm1 <- j - 1
  f[jm2] <- (cons[2] * x[jm2]^3 + cons[1] * x[jm2]) * exp(-(x[jm2]^2)/100) - 1
  f[jm1] <- 10 * (sin(x[jm2]) - x[jm1])
  f[j] <- 10 * (cos(x[jm2]) - x[j])
  sum(f*f)
  }

k <- c(1.003344481605351, -3.344481605351171e-03)
p0 <- rnorm(30)  # number of parameters should be a multiple of 3 for this example
ans.spg2 <- spg(par=p0, fn=valley.f, cons=k, method=2)  
ans.cg <- optim(par=p0, fn=valley.f, cons=k, method="CG")  
ans.lbfgs <- optim(par=p0, fn=valley.f, cons=k, method="L-BFGS-B")  

####################################################################
# Here is a statistical example illustrating log-likelihood maximization.

poissmix.loglik <- function(p,y) {
  # Log-likelihood for a binary Poisson mixture
  i <- 0:(length(y)-1)
  loglik <- y*log(p[1]*exp(-p[2])*p[2]^i/exp(lgamma(i+1)) + 
        (1 - p[1])*exp(-p[3])*p[3]^i/exp(lgamma(i+1)))
  return (sum(loglik) )
  }

# Data from Hasselblad (JASA 1969)
poissmix.dat <- data.frame(death=0:9, freq=c(162,267,271,185,111,61,27,8,3,1))
y <- poissmix.dat$freq

# Lower and upper bounds on parameters
lo <- c(0.001,0,0)  
hi <- c(0.999, Inf, Inf)

p0 <- runif(3,c(0.2,1,1),c(0.8,5,8))  # randomly generated starting values

ans.spg <- spg(par=p0, fn=poissmix.loglik, y=y, lower=lo, upper=hi, 
     control=list(maximize=TRUE))

# how to compute hessian at the MLE
  require(numDeriv)
  hess <- hessian(x=ans.spg$par, poissmix.loglik, y=y)
  se <- sqrt(-diag(solve(hess)))  # approximate standard errors


}

\keyword{multivariate}
