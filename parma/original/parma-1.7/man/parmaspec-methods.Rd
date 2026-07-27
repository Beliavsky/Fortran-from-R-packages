\name{parmaspec-methods}
\docType{methods}
\alias{parmaspec}
\alias{parmaspec-methods}
\alias{parmaspec,ANY-method}
\title{Portfolio Allocation Model Specification}
\description{
Defines the type model and method for optimization using either a scenario or 
covariance matrix.
}
\usage{
parmaspec(scenario = NULL, probability = NULL, S = NULL, Q = NULL, qB = NULL, 
		benchmark = NULL, benchmarkS = NULL, forecast = NULL, target = NULL, 
		targetType =  c("inequality", "equality"), 
		risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM", "LPMUPM"), 
		riskType = c("minrisk", "optimal", "maxreward"), riskB = NULL,
		options = list(alpha = 0.05, threshold = 999, moment = 1, lmoment=1,
		umoment=1, lthreshold = -0.01, uthreshold = 0.01), 
		LB = NULL, UB = NULL, budget = 1, leverage = NULL, 
		ineqfun = NULL, ineqgrad = NULL, eqfun = NULL, eqgrad = NULL, 
		uservars = list(), ineq.mat = NULL, ineq.LB = NULL, 
		ineq.UB = NULL, eq.mat = NULL, eqB = NULL, max.pos = NULL, 
		asset.names = NULL, ...)
}
\arguments{
\item{scenario}{ An n-by-m scenario matrix.}
\item{probability}{ An optional n-by-1 vector of scenario probabilities which
must sum to 1 (only currently used in LP problems). Default is to assign equal 
weights to each row of the scenario.}
\item{S}{ An m-by-m positive definite covariance matrix.}
\item{Q}{ A list of m-by-m positive definite matrices for QCQP type problems.}
\item{qB}{ a vector of the same length as Q denoting the upper bound on Q.}
\item{benchmark}{ A n-by-1 scenario benchmark matrix, used when scenario is not 
NULL.}
\item{benchmarkS}{ An m+1 vector consisting of the benchmark variance (1) and 
covariances (m) with the other m assets, used when S is not NULL.}
\item{forecast}{ A vector (m) of forecast values for the assets. If a benchmark
is included, this should then be the active forecast over the benchmark.}
\item{target}{ The target return required when riskType is minrisk. If a benchmark
is used, then this is the active target over the benchmark (given the active
forecasts above).}
\item{targetType}{ Whether the target should be a hard equality or inequality.}
\item{risk}{ The risk measure.}
\item{riskType}{ The type of optimization to use, with a choice of minimizing the
risk given the relevant constraints, on optimizing directly the risk to 
reward ratio using established fractional programming methods or maximizing the
reward subject to a risk upper bound and other constraints (only currently supported
by the SOCP solver for covariance matrix type problems).}
\item{riskB}{For the case that riskType is \dQuote{maxreward}, then riskB is the
upper bound for the risk constraint.}
\item{options}{ A vector of optional parameters related to the tail risk measures
CVaR, CDaR and LPM.}
\item{LB}{ The lower bounds for the asset weights. If using a benchmark, this
should be the maximum deviation below the benchmark weights.}
\item{UB}{ The upper bounds for the asset weights. If using a benchmark, this
should be the maximum deviation above the benchmark weights.}
\item{budget}{ The investment constraint. If using a benchmark, this is usually
set to zero so that weights represent the active bets on the benchmark.}
\item{leverage}{ The leverage constraint for Long/Short optimization.}
\item{ineqfun}{ A list of user inequality functions for use in an NLP type 
setup (see details).}
\item{ineqgrad}{ A list of user inequality gradients if ineqfun was provided.}
\item{eqfun}{ A list of user equality functions for use in an NLP type setup 
(see details).}
\item{eqgrad}{ A list of user equality gradients if eqfun was provided.}
\item{uservars}{ A list of any additional user required values to be used with
user defined inequality and equality functions (for NLP).}
\item{ineq.mat}{ A k-by-m inequality matrix (for LP and QP problems).}
\item{ineq.LB}{ A vector (k) of the lower bounds for the ineq.mat.}
\item{ineq.UB}{ A vector (k) of the upper bounds for the ineq.mat.}
\item{eq.mat}{ A l-by-m equality matrix (for LP and QP problems).}
\item{eqB}{ A vector (l) of equalities for the eq.mat.}
\item{max.pos}{ Cardinality Constraints. The maximum assets to include in the 
solution, effectively making this a MILP, MIQP or MINLP problem.}
\item{asset.names}{An optional character vector of asset names.}
\item{...}{ Not used. }
}
\value{
A \code{\linkS4class{parmaSpec}} object containing details of the PARMA
specification.
}
\details{
The parmaspec method is the entry point for specifying and solving portfolio 
problems in the parma package. Currently 7 measures of risk are supported, 3 
based on tail measures: Conditional Value at Risk (CVaR), Conditional Drawdown at 
Risk (CDaR) and Lower Partial Moments (LPM), and 3 based on the Lp-Norm: Mean
Absolute Deviation (L_1, MAD), Mean Variance (L_2, EV) and MiniMax 
(L_inf, Minimax). The LPMUPM measure is the ratio of lower to upper partial moments, 
a non convex measure discussed in Holthausen (1981).
Additionally, the problems may be solved based on minimization
of risk subject to a target return, else on the optimal risk-reward ratio using
fractional programming (see references), thus avoiding the estimation of the 
entire frontier. Problems are classified and solved according to whether they
can be formulated as Linear (LP), Mixed Integer LP (MILP), Quadratic (QP), 
Mixed Integer Quadratic (MIQP), Second Order Cone Programming (SOCP), 
Non-Linear (NLP), Mixed Integer NLP (MINLP) and Global NLP (GNLP). This in turn
depends on the intersection of objectives and constraints. It is possible that a
problem may be solved both as an LP and NLP (or QP and NLP), and this can be 
defined during the solver stage (\code{\link{parmasolve}}). Because all NLP models, 
make use of analytical derivatives, the results should be the same for any 
formulation chosen, and considerations such as memory usage should guide the 
choice of formulation (with some LP models being particularly expensive). Not 
all problem types are supported, but this might change subject to the availability 
of solvers in R which can deal with these specific types e.g. MINLP and MIQP. 
The parmaspec also allows the input of a benchmark so that benchmark relative 
optimization is carried out.
User defined equality and inequality functions for NLP problems need to be 
properly defined to be accepted by the model, and their analytic gradients
also provided, unless the problem is solved as a GNLP in which case a derivative 
free penalty function is used. These custom constraint functions should be 
provided in a list, and should take as arguments the vector of decision 
variables \sQuote{w}, an argument called \sQuote{optvars} which is used
by the program internally, and an argument called \sQuote{uservars} which is a 
list with optional user defined values for the constraints. The examples in the 
inst folder provide some guidance, and the user is left to his own devices to 
study the underlying workings of the program to understand how to supply these.
Finally, the NLP functions which are known to be discontinuous because of the
presence of functions such as the min and max, have been re-written to take
advantage of smooth approximations to such functions, details of which may be
founds in the vignette.
The package support for GNLP is based on a choice of the cmaes solver of the 
\pkg{cmaes} package (which is not production level) or the crs solver of the 
\pkg{nloptr} package which may be defined in solver.control option of the 
\code{\link{parmasolve}} method with named argument \sQuote{solver}. High quality
GNLP solvers are not available in R and as such support for these types of 
problems is experimental at best and your mileage will vary. The problems which 
must be solved as GNLP include the \sQuote{LPMUPM} measure, all problems with
risk type \sQuote{optimal} AND cardinality constraints (\sQuote{max.pos}), and 
all problems with custom NLP constraints without derivatives, non-convex 
inequalities or non-affine equalities.
}
\author{Alexios Galanos}
\keyword{methods}
