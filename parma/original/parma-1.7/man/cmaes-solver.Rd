\name{cmaes}
\alias{cmaes}
\alias{cmaes.control}
\title{
  The Covariance Matrix Adaptation Evolution Strategy (cmaes) Solver
}
\description{
The direct translation of the Hansen's cmaes matlab code v3.60.
}
\usage{
cmaes(pars, fun, lower = rep(0, length(pars)), upper = rep(1, length(pars)), 
insigma = 1, ctrl = cmaes.control(), ...)

cmaes.control(
options = list(StopFitness = -Inf, MaxFunEvals = Inf,
MaxIter = '1e3*(N+5)^2/sqrt(popsize)', StopFunEvals = Inf,
StopIter = Inf, TolX = '1e-11*max(insigma)', TolUpX = '1e3*max(insigma)',
TolFun = 1e-12, TolHistFun = 1e-13, StopOnStagnation = TRUE,
StopOnWarnings = TRUE, StopOnEqualFunctionValues = '2 + N/3',
DiffMaxChange = Inf, DiffMinChange = 0, WarnOnEqualFunctionValues = FALSE,
EvalParallel = FALSE, EvalInitialX = TRUE, Restarts = 0,
IncPopSize = 2, PopSize = '4 + floor(3*log(N))', ParentNumber = 'floor(popsize/2)',
RecombinationWeights = c("superlinear", "linear", "constant"), 
DiagonalOnly = '0*(1+100*N/sqrt(popsize))+(N>=1000)', 
CMA = TRUE, Seed = 'as.integer(Sys.time())', DispFinal = TRUE, 
DispModulo = 100, Warnings = FALSE), 
CMA = list(cs = '(mueff+2)/(N+mueff+3)', 
damps = '1 + 2*max(0,sqrt((mueff-1)/(N+1))-1) + cs',
ccum = '(4 + mueff/N) / (N+4 + 2*mueff/N)', ccov1 = '2 / ((N+1.3)^2+mueff)',
ccovmu = '2 * (mueff-2+1/mueff) / ((N+2)^2+mueff)', active = 0))
}
\arguments{
\item{pars}{
A numeric vector of starting parameters.}
\item{fun}{
The user function to be minimized.}
\item{lower}{
A vector the lower parameter bounds.
}
\item{upper}{
A vector with the upper parameter bounds.}
\item{insigma}{
The initial coordinate wise standard deviations for the search.}
\item{ctrl}{
A list with control parameters as returned from calling the \sQuote{cmaes.control}
function.
}
\item{...}{Additional arguments passed to the user function.}
\item{options}{The main options in the cmaes.control which may be optionally
strings which are evaluated on initialization of the solver.}
\item{CMA}{The options for the active CMA.}
}
\details{
This solver has been translated from the matlab version created by Nikolaus 
Hansen and available on his website \url{http://www.cmap.polytechnique.fr/~nikolaus.hansen/cmaes_inmatlab.html}.
There is also a \pkg{cmaes} on CRAN but this does not offer the same level of
options and flexibility that the matlab version offers. For more details on
what the options mean and generally how the cmaes solver works, consult the 
relevant website and literature.
}
\author{Alexios Galanos}
\references{
Hansen, N. 2006, The CMA Evolution Strategy: A Comparing Review, \emph{Towards 
a New Evolutionary Computation (Studies in Fuzziness and Soft Computing)}, 
\bold{192}, 75--102.\cr
}
\examples{
\dontrun{
ctrl = cmaes.control()
ctrl$options$StopOnWarnings = FALSE
ctrl$cma$active = 1
ctrl$options$TolFun = 1e-12
ctrl$options$DispModulo=100
ctrl$options$Restarts = 0
ctrl$options$MaxIter = 3000
ctrl$options$TolUpX = 5
ctrl$options$PopSize = 300
test1 = cmaes(rnorm(10), fun = parma:::fsphere, 
lower = -Inf, upper = Inf, insigma = 1, ctrl = ctrl)
test2 = cmaes(rnorm(10), fun = parma:::frosenbrock, 
lower = -Inf, upper = Inf, insigma = 1, ctrl = ctrl)

ctrl = cmaes.control()
ctrl$options$StopOnWarnings = FALSE
ctrl$cma$active = 1
ctrl$options$TolFun = 1e-12
ctrl$options$DispModulo=100
ctrl$options$Restarts = 0
ctrl$options$MaxIter = 3000
ctrl$options$PopSize = 400
test3 = cmaes(rep(1, 10), fun = parma:::frastrigin10, 
lower = -50, upper = 50, insigma = 1, ctrl = ctrl)
}}
