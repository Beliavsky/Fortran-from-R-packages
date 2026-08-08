\name{psoptim-methods}
\docType{methods}
\alias{psoptim-methods}
\alias{psoptim,ANY,ANY,ANY,ANY,ANY-method}
\alias{psoptim,test.problem,missing,missing,missing,missing-method}
\title{Methods for function psoptim (Particle Swarm Optimization)}
\description{
General implementation of particle swarm optimization usable as a direct
replacement for \code{\link{optim}}.
}
\section{Methods}{
\code{signature(par = "ANY", fn = "ANY", gr = "ANY", lower = "ANY", upper = "ANY")}:

This is the standard replacement for \code{\link{optim}} without S4
object usage.

\code{signature(par = "test.problem", fn = "missing", gr = "missing",}\cr
    \code{lower = "missing", upper = "missing")}:

  This is for running PSO on a specific test problem. Typically this is
  invoked on repetitive runs of a test problem and used to assess the
  choice of parameters for the underlying PSO algorithm. The function
  is essentially a wrapper function for \code{\link{psoptim}} but
  returns an instance of \code{\linkS4class{test.result}} containing summary
  results. 
}
\keyword{methods}
\keyword{optimize}
