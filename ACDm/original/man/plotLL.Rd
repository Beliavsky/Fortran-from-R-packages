\name{plotLL}
\alias{plotLL}

\title{
Plots the response surface of the log likelihood of a fitted model.
}
\description{
Plots the log likelihood for a fitted model against either one or two of the parameters at a time. This can help to find issues with for example poor identification of a model.
}
\usage{
plotLL(fitModel, parameter1 = 1, parameter2 = NULL, 
       param1sequence = NULL, param2sequence = NULL, startpoint = NULL, 
       length.out = NULL, returnOutput = FALSE)
}

\arguments{
  \item{fitModel}{
a fitted model of class acdFit. 
}
  \item{parameter1}{
the first parameter for the log likelihood to be plotted against. Either the index of the parameter as an integer, or the name of the parameter.
}
  \item{parameter2}{
the second parameter for the log likelihood to be plotted against. Either the index of the parameter as an integer, or the name of the parameter. If left empty, a plot with only the \code{parameter1} will be drawn.
}
  \item{param1sequence, param2sequence}{
the sequence of points from which the log likelihood is computed. If left empty, the log likelihood will be computed at 21 points spanning between MLE-3*SD and MLE+3*SD in the one dimensional case, and the 11x11 points for the same range in the two dimensional case.
}
  \item{startpoint}{
a vector of size equal to the number of parameters in the model. If this is supplied, the log likelihood will be evaluated at this point instead of the point of the MLE (for the parameters not in \code{parameter1} and \code{parameter2}).
}
  \item{length.out}{
(optional) overrides the default number of points (for each dimension) where the log likelihood is computed.
}
  \item{returnOutput}{
a logical flag. If set to \code{TRUE}, the values of the response surface will be returned. See 'value' below.
}
}

\value{
Only if \code{returnOutput = TRUE}

1. For the one dimensional case: a data.frame with the columns 'logLikelihood', and 'param1sequence' for all the values of the \code{parameter1} witch the log likelihood was evaluated at

2. For the two dimensional case: 
a list with the following items:

\item{para1}{
        a vector with the sequence of the \code{parameter1} values.
    }
\item{para2}{
        a vector with the sequence of the \code{parameter2} values.
    }
\item{z}{
        a matrix with the log likelihood values. The element at the ith row and jth column is evaluated at the ith \code{para1} value and jth \code{para2} value.
    }
    

}
\author{
Markus Belfrage
}

\examples{

\dontrun{


data(adjDurData)
#Indicates identification issues with the generelized gamma distibution:
#(Try a diffrent 'startPara' in acdFit() to get slightly a better fit)
fit <- acdFit(durations = adjDurData[1:3000, ], dist = "gengamma")
seq1 <- seq(500, 1000, 50)
seq2 <- seq(.02, 0.045, 0.001)
plotLL(fitModel = fit, parameter1 = "kappa", parameter2 = "gamma", 
       param1sequence = seq1, param2sequence = seq2)
}
}

