\name{pulmonary}
\alias{pulmonary}
\docType{data}
\title{ Change in Pulmonary Response after Exposure to Cotton Dust}
\description{
  Changes in pulmonary function of 12 workers after 6 hours of exposure to cotton dust. 
}
\usage{data(pulmonary)}
\format{
  A data frame with 12 observations on the following 3 variables.
  \describe{
    \item{\code{FVC}}{change in FVC (forced vital capacity) after 6 hours.}
    \item{\code{FEV}}{change in FEV_3 (forced expiratory volume) after 6 hours.}
    \item{\code{CC}}{change in CC (closing capacity) after 6 hours.}
  }
}
\note{
There is also a different version of this data set around. In the different version the FVC value of subject 11 is -0.01 instead of -0.10.
}

\source{
  \cite{Merchant, J. A., Halprin, G. M., Hudson, A. R. Kilburn, K. H., McKenzie, W. N., Hurst, D. J. and Bermazohnm P. (1975), 
  Responses to cotton dust, \emph{Archives of Environmental Health}, \bold{30},  222--229, Table 5. } 
  
  Reprinted with permission of the Helen Dwight Reid Educational Foundation.  
  Published by Heldref Publications, 1319 Eighteenth St., NW, Washington, DC 20036-1802. 
}
\references{
  \cite{Hettmansperger, T. P. and McKean, J. W. (1998), \emph{Robust Nonparametric Statistical Methods}, London: Arnold.  }}
  
\examples{
data(pulmonary)
plot(pulmonary)
}
\keyword{datasets}
