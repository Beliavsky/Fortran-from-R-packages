\name{getGlobalOpt}
\alias{getGlobalOpt}
\title{Return the global optimum of a given objective function. 
}
\description{Return the real value representing the global optimum
  (minimum) of a given objective function between the default bounds
  possible to obtain via the function \code{getDefaultBounds(fnName)}. 
}
\usage{
getGlobalOpt(fnName)
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
\value{Real representing the global optimum (minimum) of the objective
  function named in \code{fnName}.
}
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
\author{Katharine Mullen \email{katharine.mullen@stat.ucla.edu} }
\examples{
getGlobalOpt("Ackleys")
getGlobalOpt("AluffiPentini")
}
\keyword{optimization}