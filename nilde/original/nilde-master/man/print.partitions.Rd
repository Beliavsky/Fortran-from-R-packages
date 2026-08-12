\name{print.partitions}
\alias{print.partitions}
\alias{print.nlde}
\alias{print.knapsack}
\alias{print.subsetsum}
\alias{print.tsp_solver}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Print partitions object.}
\description{ The default print method for a \code{partitions}, \code{nlde} objects. 
}
\usage{
\method{print}{partitions}(x, ...)
\method{print}{nlde}(x, ...)
\method{print}{knapsack}(x, ...)
\method{print}{subsetsum}(x, ...)
\method{print}{tsp_solver}(x, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
\item{x, ...}{ objects of class \code{partitions}, \code{nlde}, \code{knapsack}, \code{subsetsum} as produced by \code{get.partitions()}, 
    \code{nlde()}, \code{get.knapsack()}, \code{get.subsetsum()}, \code{tsp_solver()} correspondingly.}
}
\details{ Prints the number of partitions/solutions obtained and all resulted partitions/solutions themselves. 

}


\author{
Natalya Pya Arnqvist <nat.pya@gmail.com>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%    \code{\link{S2.test}}, \code{}    
%}

%\examples{

%}











