\name{parmaSpec-class}
\docType{class}
\alias{parmaSpec-class}
\alias{show,parmaSpec-method}
\alias{parmaset<-}
\alias{parmaset<-,ANY,ANY-method}
\alias{parmaset<-,parmaSpec,vector-method}
\alias{parmaget,parmaSpec-method}
\alias{parmaget,ANY-method}
\alias{parmaget}
\title{Class \code{"parmaSpec"}}
\description{
Object returned from calling \code{\link{parmaspec}}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("parmaSpec")}.
}
\section{Slots}{
  \describe{
    \item{\code{model}:}{A list with details of the risk and optimization model.}
    \item{\code{modeldata}:}{A list with the data.}
    \item{\code{constraints}:}{A list with details on the optimization 
    constraints.}
  }
}
\section{Methods}{
  \describe{
    \item{parmasolve}{\code{signature(spec = "parmaSpec")}: 
    Solves for the optimal weights. }
    \item{show}{\code{signature(object = "parmaSpec")}: Summary method. }
    \item{parmaset<-}{\code{signature(object = "parmaSpec")}: Set a specified 
    \sQuote{arg} from the list of arguments defined in the \code{\link{parmaspec}} 
    to some \sQuote{value} and return a modified specification. }
    \item{parmaget}{\code{signature(object = "parmaSpec")}:  Get a specified 
    \sQuote{arg} from the list of arguments defined in the \code{\link{parmaspec}}.}
	 }
}

\author{Alexios Galanos}
\examples{
showClass("parmaSpec")
}
\keyword{classes}
