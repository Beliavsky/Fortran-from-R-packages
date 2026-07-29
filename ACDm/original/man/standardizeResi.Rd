\name{standardizeResi}
\alias{standardizeResi}

\title{
Residual standardization
}
\description{
Standardizes residuals from a fitted ACD model of class 'acdFit' by a probability integral transformation (taking the CDF, using the estimated distribution parameters, of the residuals) or by returning the Cox-Snell residuals.
}
\usage{
standardizeResi(fitModel, transformation = "probIntegral")
}
\arguments{
  \item{fitModel}{
a fitted ACD model of class 'acdFit'.
}
  \item{transformation}{
type of transformation done, either "probIntegral", or "cox-snell".
}
}

\details{
The probability integral transformation is done by taking the CDF of the residuals from the model estimation, using the estimated distribution parameters. Under correct specification the probability integral transformed residuals should be iid. uniform(0, 1). 

The Cox-Snell residuals is the computed by taking the integrated hazard of the residuals from the model estimation, using the estimated distribution parameters. Under correct specification the probability integral transformed residuals should be iid. unit exponentially distributed. 
}