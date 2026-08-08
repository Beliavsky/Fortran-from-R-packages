\name{psoptim}
\alias{psoptim}
\title{Particle Swarm Optimizer}
\description{
General implementation of particle swarm optimization usable as a direct
replacement for \code{\link{optim}}.
}
\usage{
psoptim(par, fn, gr = NULL, ..., lower = -1, upper = 1, control = list())
}
\arguments{
  \item{par}{Vector with length defining the dimensionality of the
    optimization problem. Providing actual values of \code{par} are not
    necessary (\code{NA} is just fine). Included primarily for
    compatibility with
    \code{\link{optim}} but if values are provided within the \code{lower}
      and \code{upper} bounds then the first particle will be
      initialized to the position provided by \code{par}.
  }
  \item{fn}{
    A function to be minimized (or maximized), with first
    argument the vector of parameters over which minimization is
    to take place.  It should return a scalar result.
  }
  \item{gr}{
    A function to return the gradient if local search is BFGS.
    If it is \code{NULL}, a
    finite-difference approximation will be used.
  }
  \item{\dots}{
    Further arguments to be passed to \code{fn} and \code{gr}.
  }
  \item{lower}{
    Lower bounds on the variables.
  }
  \item{upper}{
    Upper bounds on the variables.
  }
  \item{control}{
    A list of control parameters. See \dQuote{Details}.
  }
}
\details{
  By default this function performs minimization using a particle swarm
  algorithm, but it will maximize if \code{control$fnscale} is negative.

  The default control arguments implies that the algorithm follows the
  Standard PSO 2007 implementation by Maurice Clerc, but the code also
  provides support for PSO 2011, clamping the maximal velocity, restarting when
  all particles converge to a single area and using BFGS as
  the local search direction.

  The \code{control} argument is a list that can supply any of the
  following components:
  \describe{
  \item{trace:}{Non-negative integer. If positive, tracing information on
    the progress of the optimization is produced. Defaults to \code{0}.
  }
  \item{fnscale:}{An overall scaling to be applied to the value of \code{fn}
    and \code{gr} (if used) during optimization. If negative, turns the problem
    into a maximization problem. Optimization is performed on
    \code{fn(par)/fnscale}. Defaults to \code{1}.
  }
  \item{maxit:}{
    The maximum number of iterations. Defaults to \code{1000}.
  }
  \item{maxf:}{
    The maximum number of function evaluations (not considering any
    performed during numerical gradient computation). Defaults to \code{Inf}.
  }
  \item{abstol:}{
    The absolute convergence tolerance. The method converges once the
    best fitness obtained is less than or equal to
  \code{abstol}. Defaults to \code{-Inf}.
  }
  \item{reltol:}{
    The tolerance for restarting. Once the maximal distance between the
    best particle and all other particles is less than \code{reltol*d}
    the algorithm restarts. Defaults to \code{0} which disables the
    check for restarting.
  }
  \item{REPORT:}{
    The frequency for reports if \code{control$trace} is
    positive. Defaults to \code{10}.
  }
  \item{trace.stats:}{Logical; if \code{TRUE} statistics at every
  reporting step are collected and returned. Defaults to \code{FALSE}.
  }
  \item{s:}{
    The swarm size. Defaults to \code{floor(10+2*sqrt(length(par)))} unless \code{type} is \dQuote{SPSO2011} in which case the default is \code{40}.
  }
  \item{k:}{
    The exponent for calculating number of informants. Defaults to \code{3}.
  }
  \item{p:}{
    The average percentage of informants for each particle. A value of
    \code{1} implies that all particles are fully informed. Defaults to
    \code{1-(1-1/s)^k}.
  }
  \item{w:}{
    The exploitation constant. A vector of length \code{1} or
    \code{2}. If the length is two, the actual constant used is gradially
    changed from \code{w[1]} to \code{w[2]} as the number of iterations or
    function evaluations approach the limit provided.
    Defaults to \code{1/(2*log(2))}.
  }
  \item{c.p:}{
    The local exploration constant. Defaults to \code{.5+log(2)}.
  }
  \item{c.g:}{
    The global exploration constant. Defaults to \code{.5+log(2)}.
  }
  \item{d:}{
    The diameter of the search space. Defaults to the euclidean distance
    between \code{upper} and \code{lower}.
  }
  \item{v.max:}{
    The maximal (euclidean) length of the velocity vector. Defaults to \code{NA}
    which disables clamping of the velocity. However, if specified the
    actual clamping of the length is \code{v.max*d}.
  }
  \item{rand.order:}{
    Logical; if \code{TRUE} the particles are processed in random
    order. If \code{vectorize} is \code{TRUE} then the value of
    \code{rand.order} does not matter. Defaults to \code{TRUE}.
  }
  \item{max.restart:}{
    The maximum number of restarts. Defaults to \code{Inf}.
  }
  \item{maxit.stagnate:}{
    The maximum number of iterations without improvement.
    Defaults to \code{Inf}.
  }
  \item{vectorize:}{
    Logical; if \code{TRUE} the particles are processed in a vectorized
    manner. This reduces the overhead associated with iterating over
    each particle and may be more time efficient for cheap function
    evaluations. Defaults to \code{FALSE}.
  }
  \item{hybrid:}{
    If true, each normal PSO position update is followed by an
    L-BFGS-B search with the provided position as initial guess. This
    makes the implementation a hybrid approach. Defaults to
    \code{FALSE} which disables BFGS for the local search. Note that
    no attempt is done to control the maximal number of function
    evaluations within the local search step (this can be done
    separately through \code{hybrid.control}) but the number of
    function evaluations used by the local search method counts
    towards the limit provided by \code{maxf} AFTER the local search
    returns. To support a broader class of hybrid approaches a
    character vector can also be supplied with \dQuote{off} being
    equivalent to false, \dQuote{on} equivalent to true, and
    \dQuote{improved} implying that the local search will only be
    performed when the swarm finds an improvement.
  }
  \item{hybrid.control:}{
    List with any additional control parameters to pass on to
    \code{\link{optim}} when using L-BFGS-B for the local search.
    Defaults to \code{NULL}.
  }
  \item{type:}{Character vector which describes which reference implementation of SPSO is followed. Can take the value of \dQuote{SPSO2007} or \dQuote{SPSO2011}. Defaults to \dQuote{SPSO2007}.}
}
}
\value{
  A list, compatible with the output from \code{\link{optim}}, with components:
  \item{par}{
    The best set of parameters found.
  }
  \item{value}{
    The value of \code{fn} corresponding to \code{par}.
  }
  \item{counts}{
    A three-element vector containing the number of function
    evaluations, the number of iterations, and the number of restarts.
  }
  \item{convergence}{
    An integer code. \code{0} indicates that the algorithm terminated by
    reaching the absolute tolerance; otherwise:
    \describe{
    \item{\code{1}:}{Maximal number of function evaluations reached.}
    \item{\code{2}:}{Maximal number of iterations reached.}
    \item{\code{3}:}{Maximal number of restarts reached.}
    \item{\code{4}:}{Maximal number of iterations without improvement reached.}}
  }
  \item{message}{A descriptive message of the reason for termination.
  }

If \code{trace} is positive and \code{trace.stats} is \code{TRUE} additionally the component:
\item{stats}{A list of statistics collected at every reporting step with the following components:
\describe{\item{\code{it}}{A vector with the iteration numbers}
\item{\code{error}}{A vector with the corresponding best fitness values obtained}
\item{\code{f}}{A list with the corresponding current swarm fitness values as a vector}
\item{\code{x}}{A list with the corresponding current swarm positions as a matrix}
}}
}
\references{
  Default parameters follow:

  Clerc, M. (2011)
  \url{https://hal.archives-ouvertes.fr/hal-00764996/document}.
  Notice that the SPSO 2011 implementation does not include any
  of the bells and whistles from the implementation by M. Clerc et
  al. and effectively only differes from the SPSO 2007 implementation
  in the default swarm size, how velocities are initiated and the
  update of velocities/positions which in the SPSO 2011 implementation
  are invariant to rotation.

The gradual change of \code{w} and clamping the maximal velocity is
described in:

Parsopoulos, K.E. and Vrahatis M.N. (2002)
\emph{Recent approaches to global optimization problems
through Particle Swarm Optimization}.
Natural Computing 1: 235-306.


The restart (provided through \code{reltol}) is similar to:

Evers G.I. and Ghalia M.B.
\emph{Regrouping Particle Swarm Optimization:
A New Global Optimization Algorithm with
Improved Performance Consistency Across
Benchmarks}. \url{https://bee22.com/resources/Evers\%202009.pdf}

  The hybrid approach is similar to:

  Qin J., Yin Y. and Ban X. (2010)
  \emph{A Hybrid of Particle Swarm Optimization and Local Search for Multimodal Functions}.
  Lecture Notes in Computer Science, Volume 6145/2010, 589-596, DOI: 10.1007/978-3-642-13495-1_72
}

\seealso{
\code{\link{optim}}, \code{\linkS4class{test.problem}}.
}
\examples{
set.seed(1)
## Rastrigin function
psoptim(rep(NA,2),function(x) 20+sum(x^2-10*cos(2*pi*x)),
        lower=-5,upper=5,control=list(abstol=1e-8))

set.seed(1)
## Rastrigin function - local refinement with L-BFGS-B on improvements
psoptim(rep(NA,2),function(x) 20+sum(x^2-10*cos(2*pi*x)),
        lower=-5,upper=5,control=list(abstol=1e-8,hybrid="improved"))

## Griewank function
psoptim(rep(NA,2),function(x) sum(x*x)/4000-prod(cos(x/sqrt(1:2)))+1,
        lower=-100,upper=100,control=list(abstol=1e-2))

set.seed(1)
## Rastrigin function with reporting
o <- psoptim(rep(NA,2),function(x) 20+sum(x^2-10*cos(2*pi*x)),
             lower=-5,upper=5,control=list(abstol=1e-8,trace=1,REPORT=1,
             trace.stats=TRUE))
\dontrun{
plot(o$stats$it,o$stats$error,log="y",xlab="It",ylab="Error")
points(o$stats$it,sapply(o$stats$f,min),col="blue",pch=2)
}
}
\keyword{ optimize }
