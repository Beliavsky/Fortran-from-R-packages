\name{parmautility-methods}
\docType{methods}
\alias{parmautility}
\alias{parmautility-methods}
\alias{parmautility,ANY-method}
\title{Utility Based Optimization}
\description{
Utility based portfolio optimization using either Taylor series expansion of
utility function with moments or scenario based.
}
\usage{
parmautility(U = c("CARA", "Power"), method = c("moment", "scenario"),
		scenario = NULL, M1 = NULL, M2 =  NULL, M3 = NULL, M4 = NULL, RA = 1, 
		budget = 1, LB = rep(0, length(M1)), UB = rep(1, length(M1)))
}
\arguments{
\item{U}{ The utility function (only CARA curretly implemented).}
\item{method}{ Whether to use moment or scenario based optimization (only moment
currently implemented).}
\item{scenario}{ A n-by-m scenario matrix.}
\item{M1}{ A vector (m) of forecasts.}
\item{M2}{ An m-by-m positive definite covariance matrix.}
\item{M3}{ An m-by-m^2 third co-moment matrix.}
\item{M4}{ An m-by-m^3 fourth co-moment matrix.}
\item{RA}{Risk Aversion Coefficient for CARA.}
\item{budget}{ The investment constraint.}
\item{LB}{ The lower bounds for the asset weights (positive).}
\item{UB}{ The upper bounds for the asset weights.}
}
\value{
A \code{\linkS4class{parmaPort}} object containing details of the PARMA
optimized portfolio.
}
\details{
The function currently only implements the CARA moment based approach, but will
be expanded in the future. The moment approach can take as inputs either M1 and 
M2 (2-moment approximation), or M1, M2, M3 and M4 (4-moment approximation). Not
many models generate M3 and M4, but the \dQuote{gogarch} model with manig or 
magh distribution will.
}
\author{Alexios Galanos}
\references{
Galanos, A. and Rossi, E. and Urga, G. 2012, Independent Factor Autoregressive
Conditional Density Model \emph{submitted-TBA}\cr
}
\keyword{methods}
