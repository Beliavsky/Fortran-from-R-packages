\name{parmaPort-class}
\docType{class}
\alias{parmaPort-class}
\alias{show,parmaPort-method}
\alias{weights,parmaPort-method}
\alias{tictoc,parmaPort-method}
\alias{tictoc,ANY-method}
\alias{tictoc}
\alias{checkarbitrage,parmaPort-method}
\alias{checkarbitrage,ANY-method}
\alias{checkarbitrage}
\alias{parmarisk,parmaPort-method}
\alias{parmarisk,ANY-method}
\alias{parmarisk}
\alias{parmareward,parmaPort-method}
\alias{parmareward,ANY-method}
\alias{parmareward}
\alias{parmastatus,parmaPort-method}
\alias{parmastatus,ANY-method}
\alias{parmastatus}
\title{Class \code{"parmaPort"}}
\description{
The parma optimal portfolio class.
}
\section{Objects from the Class}{
Objects can be created by calls to \code{\link{parmasolve}.}.
}
\section{Slots}{
  \describe{
    \item{\code{solution}:}{Object of class \code{"vector"} The list with the
    optimal values.}
    \item{\code{model}:}{Object of class \code{"vector"} A list with details of 
    the risk model.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "parmaPort")}: Summary. }
    \item{weights}{\code{signature(object = "parmaPort")}: Extracts the optimal 
    weights.}
    \item{tictoc}{\code{signature(object = "parmaPort")}: Extracts the time 
    elapsed to solve the problem.}
    \item{checkarbitrage}{\code{signature(object = "parmaPort")}: Extracts the 
    arbitrage check on the scenario.}
     \item{parmarisk}{\code{signature(object = "parmaPort")}: Extracts the expected 
     risk of the optimized portfolio.}
     \item{parmareward}{\code{signature(object = "parmaPort")}: Extracts the expected 
     reward of the optimized portfolio.}
      \item{parmastatus}{\code{signature(object = "parmaPort")}: Solver termination
      code for the LP and NLP solvers.}
	 }
}
\author{Alexios Galanos}
\examples{
showClass("parmaPort")
}
\keyword{classes}
