\name{getProblemDimen}
\alias{getProblemDimen}
\title{Get the length of the parameter vector expected by a given objective function.
}
\description{This function returns an integer value that corresponds to
  the length of
  the parameter vector expected by a given objective function.  
}
\usage{
getProblemDimen(fnName)
}
\arguments{
  \item{fnName}{A character vector naming the objective function to
    call. Options are elements of
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
}
\value{An integer value
  describing the length of the parameter vector expected by
  the objective function.  
}
\details{Note that some of the above functions can handle parameter
  vectors of multiple lengths; the return value given indicates the
  parameter vector length that corresponds to a problem dimension for
  which the function \code{\link{getGlobalOpt}}
  returns the global optimum.  } 
\author{Katharine Mullen \email{katharine.mullen@stat.ucla.edu} }
\references{
  Mullen, K. M. (2014), Continuous Global Optimization in R,
  \emph{Journal of Statistical Software}, \bold{60}, 6, 1--45,
  URL \url{http://www.jstatsoft.org/v60/i06/}.
}

\examples{
getProblemDimen("Ackleys")
getProblemDimen("AluffiPentini")
}
\keyword{optimization}