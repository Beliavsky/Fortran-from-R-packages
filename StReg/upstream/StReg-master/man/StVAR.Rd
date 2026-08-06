\name{StVAR}
\alias{StVAR}

\alias{StVAR}

\docType{package}

\title{
Student's t Vector Autoregression (StVAR)
}

\description{
Maximum likelihood estimation of Student's t vector autoregression (VAR) model is the purpose of this function. It can be used to estimate the linear autoregressive function (conditional mean) and the quadratic 
autoskedastic function (conditional variance). Users can specify the model with deterministic variables such as trends and dummies in the matrix form.
}


\usage{
StVAR(Data,Trend=1,lag=1,v=1,maxiter=1000,meth="BFGS",hes="FALSE",init="na")
}


\arguments{
  \item{Data}{A data matrix with at least two columns. Cannot be empty.}

  \item{Trend}{A matrix with columns representing deterministic variables like trends and dummies. If 1 (default), model with only constant intercept is estimated. If 0, the model is estimated without the intercept term.}

  \item{lag}{A positive integer (default value is 1) as lag length.}

  \item{v}{A scalar (default value is 1) greater than or equal to 1. Degrees of freedom parameter.}

  \item{maxiter}{Maximum number of iteration. Must be an integer bigger than 10.}

  \item{meth}{One of the optimization method from \code{optim} function (default value is BFGS). See details of \code{optim} function.}

  \item{hes}{Logical (default value is FALSE). If TRUE produces estimated hessian matrix and the standard errors of estimates.}

  \item{init}{If na (default), initial values for optimization are generated from a uniform distribution. A vector of initial values can also be used
(not recommended). The length of the init vector must be equal to the number of parameters of the joint distribution.}

}

\details{For the functional form of the autoregressive function and the autoskedastic function, see Spanos (1994) and Poudyal (2012).}

\value{
\item{beta}{coefficients of the autoregressive function including the coefficients of trends in matrix form with standard errors and p-values. If some of the standard errors are 
NA's, the \code{StVAR()} function has to be run again.}
\item{var.coef}{coefficients of the autoskedastic (conditional variance) function, standard errors and p-values.}
\item{like}{maximum log likelihood value.}
\item{sigma}{contemporary variance-covariance matrix.}
\item{cvar}{\code{(v/(v+lag*l-2))*sigma*cvar} is the fitted value of the autoskedastic function where \code{l} is the rank of \code{Data}.}
\item{trends}{estimated trends in the variables.}
\item{res}{non-standardized residuals.}
\item{fitted}{fitted values of the autoregressive function.}
\item{init}{estimates of the joint distribution parameters. It can be used as new initial value \code{init} in \code{StVAR()} to improve optimization further.}
\item{hes}{the estimated hessian matrix if \code{hes=TRUE}.}
\item{S}{variance covariance matrix of the joint distribution.}
\item{R.squared}{R-squared of each equation of the model.}
\item{F.stat}{F-statistics of each equation of the model.}
\item{ad}{Anderson-Darling test for Student's t distribution.}
}


\author{
Niraj Poudyal \email{niraj.poudyal@ku.edu.np}
}


\references{
Poudyal, N. and Spanos, A. (2022), Model Validation and DSGE Modeling. Econometrics, 10 (2), 17.

Spanos, A. (1994), On Modeling Heteroskedasticity: the Student's t and Elliptical
Linear Regression Models. Econometric Theory, 10: 286-315.
}


\examples{
 ## StVAR Model#####
  ## Random number seed
  set.seed(7504)

  ## Creating trend variable.
  t <- seq(1,50,1)
  y <- x <- vector(length=50)
  y[1] <- x[1] <- 0
  # Generating data on y and x.
  ut <- rt(50,df=5)
  for(i in 2:50)
  {
  y[i] <-  0.004 + 0.0045*t[i] - 0.9*y[i-1] + .2*x[i-1] + ut[i]
  x[i] <-  0.05 - 0.005*t[i] + 0.8*x[i-1] + 1*y[i-1] + ut[i]
  }

  # The trend matrix

  Trend <- cbind(1,t)

  # Estimating the model
  stvar <- StVAR(cbind(y,x),lag=1,Trend=Trend,v=6,maxiter=2000,hes=TRUE)
  
  # Generate arbitrary dates
  dates <- seq(as.Date("2014/1/1"), as.Date("2016/1/1"), "weeks")

  ## Plotting the variable y, its estimated trend and the fitted value.
  oldpar <- par(mfcol=c(3,1))
  matplot(dates[2:length(y)],cbind(y[2:length(y)],stvar$fit[,1],stvar$trend[,1]),xlab="Months",
  type='l',lty=c(1,2,3), lwd=c(1,1,3),col=c("black","blue","black"),ylab=" ",xaxt="n")
  axis.Date(1, at = seq(as.Date("2014/1/1"), as.Date("2016/1/1"), "months"),labels=TRUE)
  legend("bottomleft",legend=c("data","trend","fitted values"),lty=c(1,2,3),lwd=c(1,1,3),
  col=c("black","blue","black"),cex=.85)
  hist(stvar$res[,1],main="",xlab="") ## Histogram of y
  matplot(dates[2:length(y)],cbind(stvar$cvar),xlab="Months",type='l',lty=2,lwd=1,
  ylab="fitted variance",xaxt="n")
  axis.Date(1, at = seq(as.Date("2014/1/1"), as.Date("2016/1/1"), "months"),labels=TRUE)
  par(oldpar)

  ## Plotting the variable x, its estimated trend and the fitted value.
  oldpar <- par(mfcol=c(3,1))
  matplot(dates[2:length(x)],cbind(x[2:length(x)],stvar$fit[,2],stvar$trend[,2]),
  xlab="Months",type='l',lty=c(1,2,3),lwd=c(1,1,3),col=c("black","blue","black"),ylab=" ",xaxt="n")
  axis.Date(1, at = seq(as.Date("2014/1/1"), as.Date("2016/1/1"), "months"),labels=TRUE)
  legend("bottomleft",legend=c("data","trend","fitted values"),lty=c(1,2,3),lwd=c(1,1,3),
  col=c("black","blue","black"),cex=.85)
  hist(stvar$res[,2],main="",xlab="") ## Histogram of x
  matplot(dates[2:length(x)],cbind(stvar$cvar),xlab="Months",type='l',lty=2,lwd=1,
  ylab="fitted variance",xaxt="n")
  axis.Date(1,at=seq(as.Date("2014/1/1"),as.Date("2016/1/1"),"months"),labels=TRUE)
  par(oldpar) 
}
