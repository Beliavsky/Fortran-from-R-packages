\name{colrec}
\alias{colrec}
\docType{data}
\title{Relative Survival Data}
\description{
  Survival of patients with colon and rectal cancer diagnosed in 1994-2000. 
}
\usage{data(colrec)}
\format{
  A data frame with 5971 observations on the following 7 variables:
 \describe{
\item{sex}{sex (1=male, 2=female).}    
\item{age}{age (in days).}
    \item{diag}{date of diagnosis (in date format).}
\item{time}{survival time (in days).}
    \item{stat}{censoring indicator (0=censoring, 1=death).}
    \item{stage}{cancer stage. Values 1-3, code \code{99} stands for unknown.}
\item{site}{cancer site. }
 }
}

\references{
Provided by Slovene Cancer Registry. The \code{age}, \code{time} and \code{diag} variables are randomly perturbed to make the identification of patients impossible.
}


\keyword{datasets}
