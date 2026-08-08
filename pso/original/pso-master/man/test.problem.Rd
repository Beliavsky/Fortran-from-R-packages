\name{test.problem}
\alias{test.problem}
\title{
Convenience constructor for the test.problem class.
}
\description{
The method enables creating of objects of class \code{"test.problem"}
for a few standard test problems.
}
\usage{
test.problem(name, n.test = 100, dim, maxf, objective, lower, upper)
}
\arguments{
  \item{name}{ The name of the test problem. Currently supports one of
    \code{"parabola"}, \code{"griewank"},
    \code{"rosenbrock"}, \code{"rastrigin"}, or \code{"ackley"}.
  }
  \item{n.test}{
    The number of tests to perform.
  }
  \item{dim}{
    Override the default dimension of the problem.
  }
  \item{maxf}{
    Override the default maximal number of function evaluations for the problem.
  }
  \item{objective}{
    Override the default objective for the function.
  }
  \item{lower}{
    Override the default lower bounds for the problem.
  }
  \item{upper}{
    Override the default upper bounds for the problem.
  }
}
\value{
  An object of class \code{"test.problem"}.
}

\seealso{
  \code{\linkS4class{test.problem}}.
}
\examples{
test.problem("rast")

test.problem("rast",dim=4,n.test=10)

test.problem("grie")
}
\keyword{datagen}
