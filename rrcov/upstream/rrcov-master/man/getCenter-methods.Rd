\name{getCenter-methods}
\docType{methods}
\title{Accessor methods to the essential slots of Cov and its subclasses}
\alias{getCenter}
\alias{getCenter-methods} 
\alias{getCov}
\alias{getCov-methods} 
\alias{getCorr}
\alias{getCorr-methods} 
\alias{getData}
\alias{getData-methods} 
\alias{getDistance}
\alias{getDistance-methods} 
\alias{getEvals}
\alias{getEvals-methods} 
\alias{getDet}
\alias{getDet-methods} 
\alias{getShape}
\alias{getShape-methods} 
\alias{getFlag}
\alias{getFlag-methods} 
\alias{getMeth}
\alias{getMeth-methods} 
\alias{isClassic}
\alias{isClassic-methods} 
\alias{getRaw}
\alias{getRaw-methods} 
\alias{names}
\alias{names-methods} 
\alias{names,Cov-method}
\alias{$-methods}
\alias{$,Cov-method}
\usage{
getCenter(obj)
getCov(obj)
getCorr(obj)
getData(obj)
getDistance(obj)
getEvals(obj)
getDet(obj)
getShape(obj)
getFlag(obj, prob=0.975)
getMeth(obj)
isClassic(obj)
getRaw(obj)
} 

\arguments{
    \item{obj}{an object of class \code{"Cov"}  or of a class derived from \code{"Cov"}}
    \item{prob}{optional argument for \code{getFlag} - probability, defaults to 0.975}
}

\description{
    Accessor methods to the slots of objects of class\code{Cov} and its subclasses
}
\section{Methods}{
\describe{

\item{obj = "Cov"}{generic functions - see \code{getCenter}, \code{getCov}, \code{getCorr}, \code{getData}, \code{getDistance}, \code{getEvals}, \code{getDet}, \code{getShape}, \code{getFlag}, \code{isClassic}}
\item{obj = "CovRobust"}{generic functions - see \code{getCenter}, \code{getCov}, \code{getCorr}, \code{getData}, \code{getDistance}, \code{getEvals}, \code{getDet}, \code{getShape}, \code{getFlag}, \code{getMeth}, \code{isClassic}}
}}
\keyword{methods}
\keyword{multivariate}
\keyword{robust}
