\name{parmafrontier-methods}
\docType{methods}
\alias{parmafrontier}
\alias{parmafrontier-methods}
\alias{parmafrontier,ANY-method}
\alias{parmafrontier,parmaSpec-method}
\title{Efficient Frontier Generator}
\description{
Solves for the portfolios on the efficient frontier given a specification
object.
}
\usage{
parmafrontier(spec, n.points = 100, miny = NULL, maxy = NULL, type = NULL, 
solver = NULL, solver.control = list(), parma.control = list(ubounds = 10000, 
mbounds = 1e+05, penalty = 10000), cluster = NULL)
}
\arguments{
\item{spec}{ A \code{\linkS4class{parmaSpec}} object.}
\item{n.points}{ The number of portfolios to solve for along the frontier.}
\item{miny}{ (Optional) Minimum return from which to calculate the frontier. If 
not provided, will be calculated by a first pass optimization.}
\item{maxy}{ (Optional) Maximum return for which to calculate the frontier. If 
not provided, will be based on the maximum forecast.}
\item{type}{ The problem type to use (the show method on \code{\linkS4class{parmaSpec}} 
will indicate the available options).}
\item{solver}{ For a GNLP type problem the choice of global solver, either
\dQuote{cmaes} or \dQuote{crs} from the \code{nloptr} package. }
\item{solver.control}{ A list with optional control parameters passed to the
\code{nloptr} or \code{\link{cmaes}} solver.}
\item{parma.control}{ Internal NLP tuning parameters, where \sQuote{ubounds} 
represents the symmetric unconstrained parameter bounds in the fractional 
problem, \sQuote{mbounds} the multiplier upper bounds in the fractional 
formulation, \sQuote{penalty} the GNLP penalty parameter.}
\item{cluster}{ A precreated cluster object from the parallel package for the
parallel evaluation of the frontier portfolios (see note). }
}
\value{
A \code{\link{matrix}} object with columns for the weights, risk, reward and
termination status of the solver (not available for the QP solver), and
n.points rows. The user should investigate the non-converged solutions as well
as remove any NA based rows (non-converged with error).
}
\details{
If using a cluster object, care should be taken since some of the LP based
problems are quite memory intensive. Though some care some gone into using
slam matrices and forcing garbage collection in intermediate setup steps of these
problems, there is no guarantee that memory will not be quickly depleted for
large problems using many cores/sockets.
}
\author{Alexios Galanos}
\keyword{methods}
