\name{griddata objects}
\alias{griddata objects}
\title{Data interpolated to a grid objects.}
\description{
These are objects of classes 

\code{griddata}

These objects store the results of interpolating a 1-D regression data set to a grid which is a power of two in length
}
\details{
The help page for \code{\link{makegrid}} and Kovac, (1997), p.81 give further details about how a \code{griddata} object is constructed. 
}
\value{
The following components must be included in a legitimate griddata object. 
\item{gridt}{a vector containing the values of the grid on the "x" axis.}
\item{gridy}{a vector containing the values of the grid on the "y" axis. This vector has to be the same length as gridt. Typically the values in (\code{gridt, gridy}) are the results of interpolating arbitrary data (\code{x,y}) onto (\code{gridt, gridy}).}
\item{G}{Codes the value of the linear interpolant matrix for the corresponding entry in \code{gridt}. The value at each point corresponds to the proportion of the original data point pointed to by \code{Gindex} that contributes to the new value at the corresponding \code{gridt} value. See Kovac, (1997), page 81 for further information.}
\item{Gindex}{Each entry in \code{Gindex} refers to one of the pairs in (\code{x,y}) which is contributing to the (\code{gridt, gridy}) interpolant. See previous help for \code{G}.} 
}
\section{GENERATION}{This class of objects is returned from the \code{\link{makegrid}} function to represent the results of interpolating a 1-D regression data set to a grid.}
\section{METHODS}{The \code{griddata} class of objects really on has one function that uses it: \code{\link{irregwd}}.}
\section{RELEASE}{
Version 3.9.6 Copyright Arne Kovac 1997
Copyright Guy Nason (help pages) 1999.  
}
\section{SEE ALSO}{
\code{\link{makegrid}}, \code{\link{irregwd}}
}
\keyword{manip}
\author{Arne Kovac}
