\name{simpleDE}
\alias{simpleDE}
\title{
  Simple Differential Evolution Algorithm
}
\description{
  Simple Differential Evolution for Minimization.
}
\usage{
simpleDE(fun, lower, upper, N = 64, nmax = 256, r = 0.4, 
            confined = TRUE, log = FALSE)
}
\arguments{
  \item{fun}{the objective function to be minimized.}
  \item{lower}{vector of lower bounds for all coordinates.}
  \item{upper}{vector of upper bounds for all coordinates.}
  \item{N}{population size.}
  \item{nmax}{bound on the number of generations.}
  \item{r}{amplification factor.}
  \item{confined}{logical; stay confined within bounds.}
  \item{log}{logical; shall a trace be printed.}
}
\details{
  Evolutionary search to minimize a function: For points in the current
  generation, children are formed by taking a linear combination of parents,
  i.e., each member of the next generation has the form
  \deqn{p_1 + r(p_2 - p_3)}
  where the \eqn{p_i} are members of the current generation and \eqn{r} is an
  amplification factor.
}
\value{
  List with the following components:
  \item{fmin}{function value at the minimum found.}
  \item{xmin}{numeric vector representing the minimum.}
  \item{nfeval}{number of function calls.}
}
\references{
  Dirk Laurie. ``A Complex Optimization". Chapter 5 In: F. Bornemann,
  D. Laurie, S. Wagon, and J. Waldvogel (Eds.). The SIAM 100-Digit Challenge.
  Society of Industrial and Applied Mathematics, 2004.
}
\author{
  HwB <hwborchers@googlemail.com>
}
\note{
  Original Mathematica version by Dirk Laurie in the SIAM textbook.
  Translated to R by Hans W Borchers.
}
\seealso{
  \code{\link{simpleEA}}, \code{DEoptim} in the `DEoptim' package.
}
\examples{
  simpleDE(fnTrefethen, lower = c(-1,-1), upper = c(1,1))
  # $fmin
  # [1] -3.306869
  # $xmin
  # [1] -0.02440308  0.21061243  # this is the true global optimum!
}
\keyword{ optimize }
