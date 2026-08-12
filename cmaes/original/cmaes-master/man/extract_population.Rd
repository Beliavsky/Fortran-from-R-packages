\name{extract_population}
\alias{extract_population}
\title{Extract the iter-th population...}
\usage{extract_population(res, iter)
}
\description{Extract the \code{iter}-th population}
\details{Return the population of the \code{iter}-th iteration of the
CMA-ES algorithm. For this to work, the populations must be saved
in the result object. This is achieved by setting
\code{diag.pop=TRUE} in the \code{control} list. Function values
are included in the result if present in the result object.}
\value{A list containing the population as the \code{par} element
and possibly the function values in \code{value} if they are
present in the result object.}
\arguments{\item{res}{A \code{cma_es} result object.}
\item{iter}{Which population to return.}
}

