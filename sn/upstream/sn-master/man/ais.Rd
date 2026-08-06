%  file sn/man/ais.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2004-2013 Adelchi Azzalini
%---------------------
\name{ais}
\alias{ais}
\docType{data}
\encoding{UTF-8}
\title{Australian Institute of Sport data}
\description{Data on 102 male and 100 female athletes collected at the 
Australian Institute of Sport, courtesy of Richard Telford and Ross Cunningham.
}
\usage{data(ais)}
\format{ 
  A data frame with 202 observations on the following 13 variables.
  \tabular{rll}{
   [,1]\tab \code{sex}\tab categorical, levels: \code{female}, \code{male}\cr
   [,2]\tab \code{sport}\tab categorical, levels: \code{B_Ball}, \code{Field},
      \code{Gym}, \code{Netball}, \code{Row}, \code{Swim}, \code{T_400m}, \cr
      \tab\tab \code{Tennis}, \code{T_Sprnt}, \code{W_Polo}\cr            
   [,3]\tab \code{RCC}\tab red cell count (numeric)\cr
   [,4]\tab \code{WCC}\tab white cell count (numeric)\cr
   [,5]\tab \code{Hc}\tab Hematocrit (numeric)\cr
   [,6]\tab \code{Hg}\tab Hemoglobin (numeric)\cr
   [,7]\tab \code{Fe}\tab  plasma ferritin concentration (numeric)\cr
   [,8]\tab \code{BMI}\tab body mass index, weight/(height)\eqn{^2}{²} 
       (numeric)\cr
   [,9]\tab \code{SSF}\tab sum of skin folds (numeric)\cr
   [,10]\tab \code{Bfat}\tab body fat percentage (numeric)\cr
   [,11]\tab \code{LBM}\tab lean body mass  (numeric)\cr
   [,12]\tab \code{Ht}\tab height, cm (numeric)\cr
   [,13]\tab \code{Wt}\tab weight, kg (numeric)\cr
  }
}
\details{The data have been made publicly available in connection with the
book by  Cook and Weisberg (1994).}

\references{
  Cook and Weisberg (1994),
  \emph{An Introduction to Regression Graphics}.                   
  John Wiley & Sons, New York.
}

\examples{
data(ais, package="sn")
pairs(ais[,c(3:4,10:13)], col=as.numeric(ais[,1]), main = "AIS data")
}
\keyword{datasets}
