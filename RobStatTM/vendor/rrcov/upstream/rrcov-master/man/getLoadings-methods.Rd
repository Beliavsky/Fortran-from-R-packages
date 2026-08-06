\name{getLoadings-methods}
\docType{methods}
\alias{getLoadings-methods}
\alias{getLoadings}
\alias{getEigenvalues,methods}
\alias{getEigenvalues}
\alias{getPrcomp,methods}
\alias{getPrcomp}
\alias{getScores,methods}
\alias{getScores}
\alias{getSdev,methods}
\alias{getSdev}
\alias{getQuan,methods}
\alias{getQuan}
\alias{getScale,methods}
\alias{getScale}
\alias{names,Pca-method}
\alias{$,Pca-method}

\title{ Accessor methods to the essential slots of Pca and its subclasses}
\description{
    Accessor methods to the slots of objects of class \code{Pca} and its subclasses
}

\arguments{
    \item{obj}{an object of class \code{"Pca"}  or of a class derived from \code{"Pca"} }
}
\section{Methods}{
\describe{
\item{obj = "Pca"}{Accessors for object of class \code{Pca}}
\item{obj = "PcaRobust"}{Accessors for object of class \code{PcaRobust}}
\item{obj = "PcaClassic"}{Accessors for object of class \code{PcaClassic}}
}}
\keyword{methods}
\keyword{multivariate}
\keyword{robust}
