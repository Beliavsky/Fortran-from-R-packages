\name{neighbourfun}
\alias{neighbourfun}
\alias{neighborfun}
\alias{neighbours}
\alias{neighbors}
\alias{neighbours-package}
\title{
  Neighbourhood Functions
}
\description{
  Create neighbourhood functions, including constraints.
}
\usage{
neighbourfun(min = 0, max = 1, kmin = NULL, kmax = NULL,
             stepsize, sum = TRUE, random = TRUE, update = FALSE,
             type = "numeric", active = TRUE, length = NULL,
             A = NULL, ...)

neighborfun (min = 0, max = 1, kmin = NULL, kmax = NULL,
             stepsize, sum = TRUE, random = TRUE, update = FALSE,
             type = "numeric", active = TRUE, length = NULL,
             A = NULL, ...)
}
\arguments{
  \item{min}{a numeric vector. A scalar is recycled to \code{length},
    if specified.
  }
  \item{max}{a numeric vector. A scalar is recycled to \code{length},
    if specified.
  }
  \item{kmin}{\code{\link{NULL}} or \code{integer}:
    the minimum number of \code{\link{TRUE}} values in logical vectors
  }
  \item{kmax}{\code{\link{NULL}} or \code{integer}:
    the maximum number of \code{\link{TRUE}} values in logical vectors
  }
  \item{stepsize}{numeric. For numeric neighbourhoods,
    the (average) stepsize. For logical neighbourhoods,
    the number of elements that are changed.
  }
  \item{sum}{logical or numeric.
    If specified and of length 1, only zero-sum
    changes will be applied to a numeric vector (i.e. the sum over all
    elements in a solution remains unchanged).
  }
  \item{random}{logical. Should the stepsize be random or fixed?
  }
  \item{active}{a vector: either the positions of elements
    that may be changed, or a logical vector. The default is a
    length-one logical vector, which means that all
    elements may be changed.
  }

  \item{update}{either \code{logical} (the default
    \code{FALSE}) or a string, specifying the type of
    updating. Currently supported is \code{"Ax"}, in which
    case the matrix \code{A} must be specified. See
    Examples.}

  \item{A}{a numeric matrix}
  \item{type}{string: either \code{"numeric"},
    \code{"logical"} or \code{"permute"}
  }
  \item{length}{
    integer: the length of a vector
  }
  \item{\dots}{
    other arguments
  }
}
\details{The function returns a closure with arguments \code{x}
  and \code{...}, which can be used for local-search
  algorithms.

  Three types of solution vectors are supported:
  \describe{

    \item{\code{numeric }}{a neighbour is created by adding
      or subtracting typically small numbers to random
      elements of a solution}%

    \item{\code{logical} }{\code{\link{TRUE}} and
      \code{\link{FALSE} } values are switched}%

    \item{\code{permute}}{elements of \code{x} are
      exchanged. Works with atomic and generic vectors (aka
      lists).}
  }%

  \code{neighborfun} is an alias for \code{neighbourfun}.

}
\value{

  A function (closure) with arguments \code{x} and
  \code{\dots}.

}
\section{Note on algorithms}{

  There are different strategies to implement constraints in
  local-search algorithms, and ultimately only experiments
  show which strategy works well for a given problem class.
  The algorithms used by \code{\link{neighbourfun}} always
  require a feasible initial solution, and then remain
  within the space of feasible solutions. See Gilli et
  al. (2019), Section 12.5, for a brief discussion.

}
\references{

  Gilli, M., Maringer, D. and Schumann, E. (2019) \emph{Numerical
  Methods and Optimization in Finance}. 2nd edition. Elsevier.\cr
  \doi{10.1016/C2017-0-01621-X}

  Schumann, E. (2025) \emph{Financial Optimisation with R
    (\pkg{NMOF} Manual)}.\cr
  \url{https://enricoschumann.net/NMOF.htm#NMOFmanual}

}
\author{
  Maintainer: Enrico Schumann <es@enricoschumann.net>
}
\seealso{%
  implementations of algorithms of the local-search family, such as
  Simulated Annealing (\code{\link[NMOF]{SAopt}} in \pkg{NMOF}) or
  Threshold Accepting (\code{\link[NMOF]{TAopt}} in \pkg{NMOF})
}
\examples{
## a LOGICAL neighbourhood
x <- logical(8)
x[1:3] <- TRUE

N <- neighbourfun(type = "logical", kmin = 3, kmax = 3)

cat(ifelse(x, "o", "."), "  | initial solution ",
    sep = "", fill = TRUE)
for (i in 1:5) {
    x <- N(x)
    cat(ifelse(x, "o", "."), sep = "", fill = TRUE)
}
## ooo.....  | initial solution
## oo....o.
## o...o.o.
## o.o.o...
## oo..o...
## oo....o.



## UPDATING a numeric neighbourhood
##   the vector is 'x' is used in the product 'Ax'
A <- array(rnorm(100*25), dim = c(100, 25))
N <- neighbourfun(type = "numeric",
                  stepsize = 0.05,
                  update = "Ax",
                  A = A)
x <- rep(1/25, 25)
attr(x, "Ax") <- A \%*\% x
for (i in 1:10)
    x <- N(x, A)

all.equal(A \%*\% x, attr(x, "Ax"))



## a PERMUTATION neighbourhood
x <- 1:5

N <- neighbourfun(type = "permute")
N(x)
## [1] 1 2 5 4 3
##         ^   ^

N <- neighbourfun(type = "permute", stepsize = 5)
N(x)

## 'x' is not restricted to integers
x <- letters[1:5]
N(x)


## a useful way to STORE/SPECIFY PARAMETERS, e.g. in config files
settings <- list(type = "numeric",
                 min = 0.0,
                 max = 0.2)
do.call(neighbourfun, settings)
}

\keyword{ optimize }
