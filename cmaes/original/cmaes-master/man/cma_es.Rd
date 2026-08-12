\name{cma_es}
\alias{cma_es}
\title{Covariance matrix adapting evolutionary strategy}
\usage{cma_es(par, fn, ..., lower, upper, control=list())
cmaES(...)
}
\description{Global optimization procedure using a covariance matrix adapting
evolutionary strategy.}
\details{\code{cma_es}: Note that arguments after \code{\dots} must be matched exactly.
By default this function performs minimization, but it will
maximize if \code{control$fnscale} is negative. It can usually be
used as a drop in replacement for \code{optim}, but do note, that
no sophisticated convergence detection is included. Therefore you
need to choose \code{maxit} appropriately.

If you set \code{vectorize==TRUE}, \code{fn} will be passed matrix
arguments during optimization. The columns correspond to the
\code{lambda} new individuals created in each iteration of the
ES. In this case \code{fn} must return a numeric vector of
\code{lambda} corresponding function values. This enables you to
do up to \code{lambda} function evaluations in parallel.

The \code{control} argument is a list that can supply any of the
following components:
\describe{
\item{\code{fnscale}}{An overall scaling to be applied to the value
of \code{fn} during optimization. If negative,
turns the problem into a maximization problem. Optimization is
performed on \code{fn(par)/fnscale}.}
\item{\code{maxit}}{The maximum number of iterations. Defaults to
\eqn{100*D^2}, where \eqn{D} is the dimension of the parameter space.}
\item{\code{stopfitness}}{Stop if function value is smaller than or
equal to \code{stopfitness}. This is the only way for the CMA-ES
to \dQuote{converge}.}
\item{keep.best}{return the best overall solution and not the best
solution in the last population. Defaults to true.}
\item{\code{sigma}}{Initial variance estimates. Can be a single
number or a vector of length \eqn{D}, where \eqn{D} is the dimension
of the parameter space.}
\item{\code{mu}}{Population size.}
\item{\code{lambda}}{Number of offspring. Must be greater than or
equal to \code{mu}.}
\item{\code{weights}}{Recombination weights}
\item{\code{damps}}{Damping for step-size}
\item{\code{cs}}{Cumulation constant for step-size}
\item{\code{ccum}}{Cumulation constant for covariance matrix}
\item{\code{vectorized}}{Is the function \code{fn} vectorized?}
\item{\code{ccov.1}}{Learning rate for rank-one update}
\item{\code{ccov.mu}}{Learning rate for rank-mu update}
\item{\code{diag.sigma}}{Save current step size \eqn{\sigma}{sigma}
in each iteration.}
\item{\code{diag.eigen}}{Save current principle components
of the covariance matrix \eqn{C}{C} in each iteration.}
\item{\code{diag.pop}}{Save current population in each iteration.}
\item{\code{diag.value}}{Save function values of the current
population in each iteration.}}

}
\value{\code{cma_es}: A list with components: \describe{
\item{par}{The best set of parameters found.}
\item{value}{The value of \code{fn} corresponding to \code{par}.}
\item{counts}{A two-element integer vector giving the number of calls
to \code{fn}. The second element is always zero for call
compatibility with \code{optim}.}
\item{convergence}{An integer code. \code{0} indicates successful
convergence. Possible error codes are \describe{
\item{\code{1}}{indicates that the iteration limit \code{maxit}
had been reached.}}}
\item{message}{Always set to \code{NULL}, provided for call
compatibility with \code{optim}.}
\item{diagnostic}{List containing diagnostic information. Possible elements
are: \describe{
\item{sigma}{Vector containing the step size \eqn{\sigma}{sigma}
for each iteration.}
\item{eigen}{\eqn{d \times niter}{d * niter} matrix containing the
principle components of the covariance matrix \eqn{C}{C}.}
\item{pop}{An \eqn{d\times\mu\times niter}{d * mu * niter} array
containing all populations. The last dimension is the iteration
and the second dimension the individual.}
\item{value}{A \eqn{niter \times \mu}{niter x mu} matrix
containing the function values of each population. The first
dimension is the iteration, the second one the individual.}}
These are only present if the respective diagnostic control variable is
set to \code{TRUE}.}}

}
\seealso{\code{\link{extract_population}}}
\references{Hansen, N. (2006). The CMA Evolution Strategy: A Comparing Review. In
J.A. Lozano, P. Larranga, I. Inza and E. Bengoetxea (eds.). Towards a
new evolutionary computation. Advances in estimation of distribution
algorithms. pp. 75-102, Springer}
\author{Olaf Mersmann \email{olafm@statistik.tu-dortmund.de} and
David Arnu \email{david.arnu@tu-dortmun.de}}
\keyword{optimize}
\arguments{\item{par}{Initial values for the parameters to be optimized over.}
\item{fn}{A function to be minimized (or maximized), with first
argument the vector of parameters over which minimization is to
take place. It should return a scalar result.}
\item{\dots}{Further arguments to be passed to \code{fn}.}
\item{lower}{Lower bounds on the variables.}
\item{upper}{Upper bounds on the variables.}
\item{control}{A list of control parameters. See \sQuote{Details}.}
}
\alias{cmaES}

