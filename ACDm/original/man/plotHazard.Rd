\name{plotHazard}

\alias{plotHazard}

\title{Hazard function plot}


\description{

Estimates and plots the hazard function from an estimatated ACD model.
}


\usage{    
plotHazard(fitModel, breaks = 20, implied = TRUE, xstop)
}

\arguments{
    
    \item{fitModel}{
        an estimated model of class acdFit. Can also be a numerical vector.
    }
    \item{breaks}{
        the number of quantiles used to estimate the hazard.
    }
    \item{implied}{
        a logical flag. If TRUE then the implied hazard function using the distribution parameter estimates will be plotted together with the nonparametric estimate of the error term hazard function.
    }
    \item{xstop}{
        where to stop plotting the implied hazard.
    }

}
    

\details{
  
    This estimator of the hazard function is based on the one used by Engle and Russell (1998). It is modified sligthly to decrease its bias and inconsistency. However, the estimator is still not fully consistent when using a fixed number of breaks (quantiles). 
          
}


\references{

    Engle, R.F and Russell, J.R. (1998) Autoregressive Conditional Duration: A New Model for Irregularly Spaced Transaction Data. \emph{Econometrica}, 66(5): 1127-1162.
    
}


\author{

    Markus Belfrage
  
}

\examples{  
\dontrun{

data(adjDurData)
fitModelWei <- acdFit(adjDurData, dist = "wei")
plotHazard(fitModelWei)
}

}


