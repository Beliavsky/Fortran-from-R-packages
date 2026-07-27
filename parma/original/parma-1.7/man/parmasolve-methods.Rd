\name{parmasolve-methods}
\docType{methods}
\alias{parmasolve}
\alias{parmasolve-methods}
\alias{parmasolve,ANY-method}
\alias{parmasolve,parmaSpec-method}
\title{Portfolio Allocation Model Solver}
\description{
Solves for the optimal weights given parmaSpec defined model.
}
\usage{
parmasolve(spec, type = NULL, solver = NULL, solver.control = list(), x0 = NULL, 
w0 = NULL, parma.control = list(ubounds = 1e4, mbounds = 1e5, 
penalty = 1e4, eqSlack = 1e-05), ...)
}
\arguments{
\item{spec}{ A \code{\linkS4class{parmaSpec}} object.}
\item{type}{ The problem type to use (the show method on \code{\linkS4class{parmaSpec}} 
will indicate the available options).}
\item{solver}{ For a GNLP type problem the choice of global solver, either
\dQuote{cmaes} or \dQuote{crs} from the \code{nloptr} package. For LP type problems,
only glpk ic currently supported.} 
\item{solver.control}{ A list with optional control parameters passed to the
\code{nloptr} or \code{\link{cmaes}} solver.}
\item{x0}{ Optional starting parameters for the NLP type problems. This 
may be of size greater than the number of assets depending on whether it takes
extra optimization parameters (CVaR has VaR in position 1, whilst the fractional
formulation has the multiplier in the last position).}
\item{w0}{ Optional starting parameters for only the asset weights.}
\item{parma.control}{ Internal NLP tuning parameters, where \sQuote{ubounds} represents
the symmetric unconstrained parameter bounds in the fractional problem, \sQuote{mbounds}
the multiplier upper bounds in the fractional formulation, \sQuote{penalty} the GNLP
penalty parameter and \sQuote{eqSlack} the slack value to use for converting equalities
to inequalities in the SOCP formulation.}
\item{...}{ Currently only the \sQuote{verbose} argument (logical) which is used 
and passed to the GLPK solver (similar to \sQuote{trace}). }
}
\value{
A \code{\linkS4class{parmaPort}} object containing details of the PARMA
optimized portfolio.
}
\details{
In most of the cases, the intersection of objective and constraints will define
whether a problem is LP, MILP, QP, MIQP, QCQP, NLP, MINLP or GNLP. However, 
there are cases when the problem can be solved by more than one type of solver 
so that the \sQuote{type} option allows the user to choose between the options. 
For the GNLP type, the \sQuote{solver} allows for a choice of \dQuote{cmaes} or 
\dQuote{crs}, while the LP type there is a choice of \dQuote{glpk} or \dQuote{symphony}.
Future development will likely expand on these choices.
}
\author{Alexios Galanos}
\keyword{methods}
