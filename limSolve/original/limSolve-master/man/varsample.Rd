\name{varsample}
\alias{varsample}
\title{
  Samples the probability density function of variables of linear inverse problems.
}
\description{
  Uses random samples of an under- or overdetermined linear problem to
  estimate the distribution of equations

  Based on a random sample of x (e.g. produced with \code{\link{xsample}}),
  produces the corresponding set of "variables" consisting of linear
  equations in the unknowns.

  \deqn{Var = EqA.x+EqB}
}

\usage{
varsample (X, EqA, EqB=NULL)
}

\arguments{
  \item{X }{matrix whose rows contain the sampled values of the unknowns
    \code{x} in \eqn{EqA*x-EqB}.
  }
  \item{EqA }{numeric matrix containing the coefficients that define
    the variables.
  }
  \item{EqB }{numeric vector containing the right-hand side of the variable
    equation.
  }
}
\value{
  a matrix whose rows contain the sampled values of the variables.
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\examples{
# The probability distribution of vertebrate and invertebrate
# food in the diet of Mink
# food items of Mink are (in that order):

# fish mussels crabs shrimp rodents amphipods ducks
# V    I       I     I      V       I         V 
# V= vertebrate, I = invertebrate

# In matrix form:
VarA <- matrix(ncol = 7, byrow = TRUE, data = c(
        0, 1, 1, 1, 0, 1, 0,   # invertebrates  
        1, 0, 0, 0, 1, 0, 1))  # vertebrates

# first sample the Minkdiet problem
E <- rbind(Minkdiet$Prey, rep(1, 7))
F <- c(Minkdiet$Mink, 1)
X <- xsample(E = E, F = F, G = diag(7), H = rep(0, 7), iter = 1000)$X

#then determine Diet Composition in terms of vertebrate and invertebrate food
DC <- varsample(X = X, EqA = VarA)
hist(DC[,1], freq = FALSE, xlab = "fraction",
     main = "invertebrate food in Mink diet", col = "lightblue")
}

\seealso{
  \code{\link{Minkdiet}}, for a description of the Mink diet example.

  \code{\link{varranges}}, to estimate ranges of inverse variables.

  \code{\link{xsample}}, to randomly sample the lsei problem.
}
\keyword{algebra}
\keyword{array}

