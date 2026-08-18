\name{summarizeDF}
\alias{summarizeDF}
\title{Summarize a data.rrame}

\description{
This function summarizes each variable that is in a data.frame.
It can be used, e.g., in an R script to write summary information 
about a data.frame into a text file that is in a version control system
so that one can see in the version control system 
whether one or more variables in the data frame have changed.
}

\usage{
summarizeDF( dat, printValues = TRUE, maxLevel = 20, file = NULL, ... )
}

\arguments{
   \item{dat}{a data.frame.}
   \item{printValues}{logical. If \code{FALSE} only MD5 checksums are returned,
      which could be desirable if the data frame contains confidential data
      that should not be included in the output.}
   \item{maxLevel}{integer. If the number of unique values in a variable 
      is less than or equal to the number specified in this argument
      (and argument \code{printValues} is \code{TRUE}),
      a frequency table is included in the output.}
   \item{file}{a character string or a writable connection naming 
      the file to write to.}
   \item{...}{further arguments forwarded to \code{sink()}
      if argument \code{file} is not \code{NULL}.}
}

\author{Arne Henningsen}
\keyword{methods}
