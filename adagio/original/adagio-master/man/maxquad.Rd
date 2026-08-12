\name{maxquad}
\alias{maxquad}
\title{
  The MAXQUAD Test Function
}
\description{
  Lemarechal's MAXQUAD optimization test function.
}
\usage{
maxquad(n, m)
}
\arguments{
  \item{n}{number of variables of the generated test function.}
  \item{m}{number of functions to compete for the maximum.}
}
\details{
  MAXQUAD actually is a family of minimax functions, parametrized by the
  number \code{n} of variables and the number \code{m} of functions whose
  maximum it is.
}
\value{
  Returns a list with components \code{fn} the generated test function of
  \code{n} variables, and \code{gr} the corresponding (analytical) gradient 
  function.
}
\references{
  Kuntsevich, A., and F. Kappel (1997). SolvOpt -- The Solver for Local 
  Nonlinear Optimization Problems.
  Manual Version 1.1, Institute of Mathematics, University of Graz.

  Lemarechal, C., and R. Mifflin, Eds. (1978). Nonsmooth Optimization.
  Pergamon Press, Oxford.

  Shor, N. Z. (1985). Minimization Methods for Non-differentiable Functions.
  Series in Computational Mathematics, Springer-Verlag, Berlin.
}
\examples{
# Test function of 5 variables, defined as maximum of 5 smooth functions
maxq <- maxquad(5, 5)
fnMaxquad <- maxq$fn
grMaxquad <- maxq$gr
# shor
}

\keyword{ testfunctions }
