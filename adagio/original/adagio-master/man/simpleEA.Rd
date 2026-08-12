\name{simpleEA}
\alias{simpleEA}
\title{
  Simple Evolutionary Algorithm
}
\description{
  Simple Evolutionary Algorithm for Minimization.
}
\usage{
simpleEA(fn, lower, upper, N = 100, ..., con = 0.1, new = 0.05,
         tol = 1e-10, eps = 1e-07, scl = 1/2, confined = FALSE, log = FALSE)
}
\arguments{
  \item{fn}{the objective function to be minimized.}
  \item{lower}{vector of lower bounds for all coordinates.}
  \item{upper}{vector of upper bounds for all coordinates.}
  \item{N}{number of children per parent.}
  \item{...}{additional parameters to be passed to the function.}
  \item{con}{percentage of individuals concentrating to the best parents.}
  \item{new}{percentage of new individuals not focussed on existing parents.}
  \item{tol}{tolerance; if in the last three loops no better individuals
             were found up to this tolerance, stop.}
  \item{eps}{grid size bound to be reached.}
  \item{scl}{scaling factor for shrinking the grid.}
  \item{confined}{logical; shall the set of individuals be strictly
                  respect the boundary? Default: FALSE.}
  \item{log}{logical, should best solution found be printed per step.}
}
\details{
  Evolutionary search to minimize a function: For each point in the current
  generation, \emph{n} random points are introduced and the \emph{n} best
  results of each generation (and its parents) are used to form the next
  generation.

  The scale shrinks the generation of new points as the algorithm proceeds.
  It is possible for some children to lie outside the given rectangle, and
  therefore the final result may lie outside the unit rectangle well.
  (TO DO: Make this an option.)
}
\value{
  List with the following components:
  \item{par}{numeric vector representing the minimum found.}
  \item{val}{function value at the minimum found.}
  \item{fun.calls}{number of function calls made.}
  \item{rel.scl}{last scaling factor indicating grid size in last step.}
  \item{rel.tol}{relative tolerance within the last three minima found.}
}
\references{
  Stan Wagon. ``Think Globally, Act Locally". Chapter 4 In: F. Bornemann,
  D. Laurie, S. Wagon, and J. Waldvogel (Eds.). The SIAM 100-Digit Challenge.
  Society of Industrial and Applied Mathematics, 2004.
}
\author{
  HwB <hwborchers@googlemail.com>
}
\note{
  Original Mathematica Version by Stan Wagon in the SIAM textbook.
  Translated to R by Hans W Borchers.
}
\seealso{
  \code{DEoptim} in the `DEoptim' package.
}
\examples{
  simpleEA(fnTrefethen, lower=c(-1,-1), upper=c(1,1), log=FALSE)
  # $par
  # [1] -0.02440310  0.21061243  # this is the true global optimum!
  # $val
  # [1] -3.306869
}
\keyword{ optimize }
