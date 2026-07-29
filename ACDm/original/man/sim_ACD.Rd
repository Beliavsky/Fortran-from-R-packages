\name{sim_ACD}
\alias{sim_ACD}

\title{
ACD simulation
}
\description{
Simulates a sample from a specified ACD \code{model} and error term distribution \code{dist}. The error terms can also be sampled from residuals.
}
\usage{    
sim_ACD(N = 1000, model = "ACD", dist = "exponential", param = NULL, order = NULL,
    Nburn = 50, startX = c(1), startMu = c(1), errors = NULL, sampleErrors = TRUE, 
    roundToSec = FALSE, rm0 = FALSE, bp = NULL)
}

\arguments{
  \item{N}{
sample size
}
  \item{model}{
the class of conditional mean duration specification. One of \code{"ACD"}, \code{"LACD1"}, \code{"LACD2"}, \code{"AMACD"},\code{"ABACD"}, \code{"SNIACD"} or \code{"LSNIACD"}. Can also be an object of class \code{acdFit}, as returned by \code{acdFit()}.
}
  \item{dist}{
the distribution of the error terms (only if \code{errors} are left out). Must be one of \code{"exponential"}, \code{"weibull"}, \code{"burr"}, \code{"gengamma"} or \code{"genf"}.
}
  \item{param}{
a vector of the parameters of the DGP (data generating process).
}
  \item{order}{
a vector describing the order of the conditional mean duration specification, e.g. \code{order = c(1,1)} for an ACD(1,1) model.
}
  \item{Nburn}{
the number of burned observations. Used to lower the effect of the start values of the simulated series. 
}
  \item{startX}{
a vector of values to start the simulation from.
}
  \item{startMu}{
a vector of conditional mean values to start the simulation from.
}
  \item{errors}{
a vector of error terms. If provided and \code{sampleErrors = TRUE} the errors will be sampled from this vector (with replacement). If instead \code{sampleErrors = FALSE} the error terms will be matched by the \code{errors} vector non stochastic (must then be of the same length as \code{N + Nburn})
}
  \item{sampleErrors}{
logical flag, see \code{errors} above. Default is \code{TRUE}.
}
  \item{roundToSec}{
if \code{TRUE} the simulated sample will be discretized with 1 second(unit) precision.
}
  \item{rm0}{
if \code{TRUE} zero durations will be removed. Will the result in a smaller sample than \code{N}.
}
  \item{bp}{
a vector of breakpoints used for \code{model} = "SNIACD", "LSNIACD", "TACD", "TAMACD". 
}

}
\value{
a numerical vector of simulated ACD durations
}
\author{
Markus Belfrage
}
\examples{

x <- sim_ACD() #simulates 1000 observations from an ACD(1,1) with exp. errors as default
acdFit(x)


}