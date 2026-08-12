\name{Testfunctions}
\alias{fnRosenbrock}
\alias{grRosenbrock}
\alias{fnRastrigin}
\alias{grRastrigin}
\alias{fnNesterov}
\alias{grNesterov}
\alias{fnNesterov1}
\alias{fnNesterov2}
\alias{fnHald}
\alias{grHald}
\alias{fnShor}
\alias{grShor}
\title{
  Optimization Test Functions
}
\description{
  Simple and often used test function defined in higher dimensions and with
  analytical gradients, especially suited for performance tests. Analytical 
  gradients, where existing, are provided with the \code{gr} prefix.
  The dimension is determined by the length of the input vector.
}
\usage{
fnRosenbrock(x)
grRosenbrock(x)
fnRastrigin(x)
grRastrigin(x)
fnNesterov(x)
grNesterov(x)
fnNesterov1(x)
fnNesterov2(x)
fnHald(x)
grHald(x)
fnShor(x)
grShor(x)
}
\arguments{
  \item{x}{numeric vector of a certain length.}
}
\details{
\bold{Rosenbrock} -- Rosenbrock's famous valley function from 1960. It can
also be regarded as a least-squares problem:
\deqn{\sum_{i=1}^{n-1} (1-x_i)^2 + 100 (x_{i+1}-x_i^2)^2}
\tabular{ll}{
  No. of Vars.: \tab n >= 2 \cr
  Bounds: \tab -5.12 <= xi <= 5.12 \cr
  Local minima: \tab at f(-1, 1, \ldots, 1) for n >= 4 \cr
  Minimum: \tab 0.0 \cr
  Solution: \tab xi = 1, i = 1:n \cr
}

\bold{Nesterov} -- Nesterov's smooth adaptation of Rosenbrock, based on the 
idea of Chebyshev polynomials. This function is even more difficult to 
optimize than Rosenbrock's:
\deqn{(1 - x_1)^2 / 4 + \sum_{i=1}^{n-1} (1 + x_{i+1} - 2 x_i^2)^2}
Two nonsmooth Nesterov functions are available: Nesterov2 and Nesterov1,
defined as
\deqn{(1 - x_1)^2 / 4 + \sum_{i=1}^{n-1} |1 + x_{i+1} - 2 x_i^2|}
\deqn{|1 - x_1| / 4 + \sum_{i=1}^{n-1} (|1 + x_{i+1} - 2 |x_i||}

\tabular{ll}{
  No. of Vars.: \tab n >= 2 \cr
  Bounds: \tab -5.12 <= xi <= 5.12 \cr
  Local minima: ?\cr
  Minimum: \tab 0.0 \cr
  Solution: \tab xi = 1, i = 1:n \cr
}

\bold{Nesterov1} and \bold{Nesterov2} -- Simlar to \code{Nesterov}, except the terms added are taken with absolute value, which makes this function nonsmooth and painful for gradient-based optimization routines; no gradient provided.\cr (Nesterov2 uses absolute instead of quadratic terms.)

\bold{Rastrigin} -- Rastrigin's function is a famous, non-convex example from 1989 for global optimization. It is a typical example of a multimodal function with many local minima:
  \deqn{10 n + \sum_1^n (x_i^2 - 10 \cos(2 \pi x_i))}
\tabular{ll}{
  No. of Vars.: \tab n >= 2 \cr
  Bounds: \tab -5.12 <= xi <= 5.12 \cr
  Local minima: \tab many \cr
  Minimum: \tab 0.0 \cr
  Solution: \tab xi = 0, i = 1:n \cr
}

\bold{Hald} -- Hald's function is a typical example of a non-smooth test
function, from Hald and Madsen in 1981.
\deqn{\max_{1 \le i \le n} | \frac{x_1 + x_2 t_i}{1 + x_3 t_i + x_4 t_i^2 + x_5 t_i^3} - \exp(t_i)|}
where \eqn{n = 21} and \eqn{t_i = -1 + (i - 1)/10} for \eqn{1 \le i \le 21}.
\tabular{ll}{
  No. of Vars.: \tab n =5 \cr
  Bounds: \tab -1 <= xi <= 1 \cr
  Local minima: \tab ? \cr
  Minimum: \tab 0.0001223713 \cr
  Solution: \tab (0.99987763,  0.25358844, -0.74660757,  0.24520150, -0.03749029) \cr
}

\bold{Shor} -- Shor's function is another typical example of a non-smooth test
function, a benchmark for Shor's R-algorithm.
}
\value{
  Returns the values of the test function resp. its gradient at that point.
  If an analytical gradient is not available, a function computing the gradient 
  numerically will be provided.
}
\references{
  Search the Internet.
}
\examples{
x <- runif(5)
fnHald(x); grHald(x)

# Compare analytical and numerical gradient
shor_gr <- function(x) adagio:::ns.grad(fnShor, x)    # internal gradient
grShor(x); shor_gr(x) 
}

\keyword{ testfunctions }