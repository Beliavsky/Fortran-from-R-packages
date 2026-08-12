\name{sumt}
\alias{sumt}
\title{Sequential Unconstrained Minimization Technique}
\description{
  Solve constrained optimization problems via the Sequential
  Unconstrained Minimization Technique (\acronym{SUMT}).
}
\usage{
sumt(x0, L, P, grad_L = NULL, grad_P = NULL, method = NULL,
     eps = NULL, q = NULL, verbose = NULL, control = list())
}
\arguments{
  \item{x0}{a list of starting values, or a single starting value.}
  \item{L}{a function to minimize.}
  \item{P}{a non-negative penalty function such that \eqn{P(x)} is zero
    iff the constraints are satisfied.}
  \item{grad_L}{a function giving the gradient of \code{L}, or
    \code{NULL} (default).}
  \item{grad_P}{a function giving the gradient of \code{P}, or
    \code{NULL} (default).}
  \item{method}{a character string, or \code{NULL}.  If not given,
    \code{"CG"} is used.  If equal to \code{"nlm"}, minimization is
    carried out using \code{\link[stats]{nlm}}.  Otherwise,
    \code{\link[stats]{optim}} is used with \code{method} as the given
    method.}
  \item{eps}{the absolute convergence tolerance.  The algorithm stops if
    the (maximum) distance between successive \code{x} values is
    less than \code{eps}.

    Defaults to \code{sqrt(.Machine$double.eps)}.}
  \item{q}{a double greater than one controlling the growth of the
    \eqn{\rho_k} as described in \bold{Details}.
    
    Defaults to 10.}
  \item{verbose}{a logical indicating whether to provide some output on
    minimization progress.
    
    Defaults to \code{getOption("verbose")}.}
  \item{control}{a list of control parameters to be passed to the
    minimization routine in case \code{optim} is used.}
}
\details{
  The Sequential Unconstrained Minimization Technique is a heuristic for
  constrained optimization.  To minimize a function \eqn{L} subject to
  constraints, one employs a non-negative function \eqn{P} penalizing
  violations of the constraints, such that \eqn{P(x)} is zero iff \eqn{x}
  satisfies the constraints.  One iteratively minimizes \eqn{L(x) +
    \rho_k P(x)}, where the \eqn{\rho} values are increased according to
  the rule \eqn{\rho_{k+1} = q \rho_k} for some constant \eqn{q > 1},
  until convergence is obtained in the sense that the Euclidean distance
  between successive solutions \eqn{x_k} and \eqn{x_{k+1}} is small
  enough.  Note that the \dQuote{solution} \eqn{x} obtained does not
  necessarily satisfy the constraints, i.e., has zero \eqn{P(x)}.  Note
  also that there is no guarantee that global (approximately)
  constrained optima are found.  Standard practice would recommend to
  use the best solution found in \dQuote{sufficiently many} replications
  of the algorithm.

  The unconstrained minimizations are carried out by either
  \code{\link[stats]{optim}} or \code{\link[stats]{nlm}}, using analytic
  gradients if both \code{grad_L} and \code{grad_P} are given, and
  numeric ones otherwise.

  If more than one starting value is given, the solution with the
  minimal augmented criterion function value is returned.
}
\value{
  A list inheriting from class \code{"sumt"}, with components \code{x},
  \code{L}, \code{P}, and \code{rho} giving the solution obtained, the
  value of the criterion and penalty function at \code{x}, and the final
  \eqn{\rho} value used in the augmented criterion function. 
}

\references{
  \bibshow{Fiacco+McCormick:1968}
}
\keyword{optimize}
