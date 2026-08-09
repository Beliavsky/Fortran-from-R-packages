\name{goTest}
\alias{goTest}
\title{Call an objective function
}
\description{This function is used to access one of 50 objective
  functions that may be useful for benchmarking performance of global
  optimization algorithms. 
}
\usage{
goTest(par, fnName, checkDim = TRUE)
}
\arguments{
  \item{par}{Vector of reals representing the parameter vector at which
    to evaluate the objective function.
  }
  \item{fnName}{A character vector representing the name of the
    objective function to use.  Options are elements of 
 \preformatted{c("Ackleys", "AluffiPentini", "BeckerLago",
      "Bohachevsky1", "Bohachevsky2", "Branin",
      "Camel3", "Camel6", "CosMix2", "CosMix4",
      "DekkersAarts", "Easom", "EMichalewicz",
      "Expo", "GoldPrice", "Griewank", "Gulf",
      "Hartman3", "Hartman6", "Hosaki", "Kowalik",
      "LM1", "LM2n10", "LM2n5", "McCormic",
      "MeyerRoth", "MieleCantrell", "Modlangerman",
      "ModRosenbrock", "MultiGauss", "Neumaier2",
      "Neumaier3", "Paviani", "Periodic",
      "PowellQ", "PriceTransistor", "Rastrigin",
      "Rosenbrock", "Salomon", "Schaffer1",
      "Schaffer2", "Schubert", "Schwefel",
      "Shekel10", "Shekel5", "Shekel7",
      "Shekelfox5", "Wood", "Zeldasine10",
      "Zeldasine20").}
    Unique partial matches to these elements
    are accepted. 
  }   
\item{checkDim}{A boolean value that is \code{TRUE} if the length of the
  parameter vector should be checked to make sure it matches the
  expected value.  If set to \code{FALSE}, the function runs slightly faster.
}}
\value{A real that representes the value of the objective function
  \code{fnName} at the parameter vector \code{par}. 
}
\author{Katharine Mullen \email{katharine.mullen@stat.ucla.edu} }
\references{
  Montaz Ali, M., Khompatraporn, C. and Zabinsky, Z. B. (2005),
  A Numerical Evaluation of Several Stochastic Algorithms on Selected
  Continuous Global Optimization Test Problems, \emph{Journal
    of Global Optimization}, \bold{31}, 4, 635--672.
  
  Mullen, K. M. (2014), Continuous Global Optimization in R,
  \emph{Journal of Statistical Software}, \bold{60}, 6, 1--45,
  URL \url{http://www.jstatsoft.org/v60/i06/}.

  \url{http://www.gamsworld.org/performance/selconglobal/selcongloballib.htm}
  for the \code{C} source. 
}
\details{
  Note that the factor of 0.02 in Ackley's function has been changed to
  the value commonly found in the literature, 0.2.  Also, Storn's
  Tchebychev Problem in 9 and 17 dimensions was not included, since the
  global minima of the imlementation of these functions does not
  correspond to the global minima reported in the above paper.
}

\examples{
goTest(fnName="Ackleys", par=rnorm(10))
goTest(fnName="AluffiPentini", par=c(1,2))
goTest(fnName="AluffiPentini",
par=rep(1,getProblemDimen("AluffiPentini")))

## use in an optimization via 'optim'
optim(par=c(1,2), fn=goTest, fnName="AluffiPentini")
}
\keyword{optimization}