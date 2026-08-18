\name{hot100}
\docType{data}
\alias{hot100}
\title{Hot 100 hits}
\description{
  This data corresponds to the 965 "number 1" hits on the Hot 100 chart over the
period January 1955 to December 2003. For a recording that reaches the
number one spot, \emph{Weeks} measures the number of
weeks that it stays at number one. The covariates are: \emph{Elvis} = 1 if the
recording was by Elvis Presley, = 0 otherwise; \emph{Beatles} = 1 if the
recording was by the Beatles, = 0 otherwise; \emph{Group} = 1 if the
recording was by a band, = 0 otherwise;
\emph{Female} = 1 if the artist was a solo female, = 0
otherwise; \emph{Male} = 1 if the artist was a solo male, = 0
otherwise; \emph{Inst} = 1 if the recording was purely
instrumental, = 0 otherwise; and \emph{NonCon} = 1 if
the recording topped the charts in nonconsecutive
weeks, = 0 otherwise.  
}
\usage{hot100}
\format{A data frame with 965 rows and 9 columns.}
\source{http://web.uvic.ca/~dgiles/downloads/data/hot100.xls}
\references{
Giles, D. E. (2006) Superstardom in the US popular music industry revisited. Economics Letters, 92(1):68–74.
Giles, D. E. (2007) Modeling inflated count data. In Y. Berbers and W. Zwaenepoel, editors, Proceedings of
the MODSIM 2007 International Congress on Modelling and Simulation, pages 919–925. L. Oxley and
D. Kulasiri, Eds., Modelling and Simulation Society of Australia and New Zealand
}
\keyword{datasets}