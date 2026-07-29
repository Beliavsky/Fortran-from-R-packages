\name{ecmAsyFit}
\alias{ecmAsyFit}
\title{Fitting Asymmetric Error Correction Model}
\description{Estimate an asymmetric error correction model (ECM) for two time series.}
\usage{
  ecmAsyFit(y, x, lag = 1, split = TRUE, 
    model = c("linear", "tar", "mtar"), thresh = 0)
}
\arguments{
  \item{y}{dependent or left-side variable for the long-run regression.}
  \item{x}{independent or right-side variable for the long-run regression.}
  \item{lag}{number of lags for variables on the right side.}
  \item{split}{a logical value (default of TRUE) of whether the right-hand variables should be split into positive and negative parts.}
  \item{model}{a choice of three models: linear, tar , or mtar cointegration.}
  \item{thresh}{a threshold value; this is only required when the model is specified as 'tar' or 'mtar.'}
}
\details{
There are two specficiations of an asymmetric ECM. The first one is how to calculate the error correction terms. One way is through linear two-step Engle Granger approach, as specificied by \code{model="linear"}. The other two ways are threshold cointegration by either 'tar' or 'mtar' with a threshold value. The second specification is related to the possible asymmetric price transmission in the lagged price variables, as specified in \code{split = TRUE}. Note that the linear cointegration specification is a special case of the threshold cointegration. A model with \code{model="linear"} is the same as a model with \code{model="tar", thresh = 0}.   
}
\value{
Return a list object of class "ecm" and "ecmAsyFit" with the following components: 
  \item{y}{dependend variable} 
  \item{x}{independent variable}
  \item{lag}{number of lags}    
  \item{split}{logical value of whether the right-hand variables are split}
  \item{model}{model choice}
  \item{IndVar}{data frame of the right-hand variables used in the ECM}
  \item{name.x}{name of the independent variable}
  \item{name.y}{name of the dependent variable}
  \item{ecm.y}{ECM regression for the dependent variable }
  \item{ecm.x}{ECM regression for the independent variable}
  \item{data}{all the data combined for the ECM}
  \item{thresh}{thresh value for TAR and MTAR model}       
}

\section{Methods}{
  Two methods are defined as follows:
  \describe{
    \item{\code{print}:}{ showing the key outputs.}
    \item{\code{summary}:}{ summarizing thekey outputs.}
  }
}

\references{
Enders, W., and C.W.J. Granger. 1998. Unit-root tests and asymmetric adjustment with an example using the term structure of interest rates. Journal of Business & Economic Statistics 16(3):304-311.
}
\author{Changyou Sun (\email{edwinsun258@gmail.com})}

\seealso{\code{\link{print.ecm}}; \code{\link{summary.ecm}}; \code{\link{ecmDiag}}; and \code{\link{ecmAsyTest}}.}
\examples{
# see example at daVich
}
\keyword{regression} 