\name{StDLM}
\alias{StDLM}

\alias{StDLM}

\docType{package}

\title{
Student's t Dynamic Linear Regression Model (StDLM)
}

\description{
Maximum likelihood estimation of Student's t dynamic linear regression model is the purpose of this function. It can be used to estimate the dynamic linear autoregressive function (conditional mean) and the quadratic 
autoskedastic function (conditional variance). Users can specify the model with deterministic variables such as trends and dummies in the matrix form.
}


\usage{
StDLM(y,X,Trend=1,lag=1,v=1,maxiter=1000,meth="BFGS",hes="FALSE",init="na")
}


\arguments{
  \item{y}{A vector representing dependent variable. Cannot be empty.}

  \item{X}{A data matrix whose columns represent exogenous variables. Cannot be empty.}

  \item{Trend}{A matrix with columns representing deterministic variables like trends and dummies. If 1 (default), model with only constant intercept is estimated. If 0, the model is estimated without an intercept term.}

  \item{lag}{A non-negative integer (default value is 1) as lag length. If lag=0, a static model is estimated.}

  \item{v}{A scalar (default value is 1) greater than or equal to 1. Degrees of freedom parameter.}

  \item{maxiter}{Maximum number of iteration. Must be an integer bigger than 10.}

  \item{meth}{One of the optimization method from \code{optim} function (default value is BFGS). See details of \code{optim} function.}

  \item{hes}{Logical (default value is FALSE). If TRUE produces estimated hessian matrix and the standard errors of estimates.}

  \item{init}{If na (default), initial values for optimization are generated from a uniform distribution. A vector of initial values can also be used
(not recommended). The length of the init vector must be equal to the number of parameters of the joint distribution.}

}


\details{For the functional form of the autoregressive function and the autoskedastic function, see Spanos (1994) and Poudyal (2012).}

\value{
\item{beta}{coefficients of the dynamic linear regression model including the coefficients of trends in matrix form with standard errors and p-values. If some of the standard errors are 
NA's, the \code{StVAR()} function has to be run again.}
\item{var.coef}{coefficients of the autoskedastic function, standard errors and p-values if if \code{hes=TRUE}}
\item{like}{maximum log likelihood value.}
\item{sigma}{contemporary variance covariance matrix.}
\item{cvar}{\code{(v/(v+lag*l-2))*sigma*cvar} is the fitted value of the autoskedastic function where \code{l} is the rank of \code{Data}}
\item{trend}{estimated trend in the variable y.}
\item{res}{non-standardized residuals}
\item{fitted}{fitted values of the autoregressive function.}
\item{init}{estimates of the joint distribution parameters. It can be used as new initial value \code{init} in \code{StVAR()} to improve optimization further.}
\item{S}{variance covariance matrix of the joint distribution.}
\item{R.squared}{R-squared of the model.}
\item{F.stat}{F-statistics of the model.}
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
## StDLM Model#####
## Random number seed
set.seed(7504)

## Creating trend variable.
t <- seq(1,100,1)
ut <- rt(100,df=5)
# Generating data on y, x and z. 
y <-  0.004 + 0.0045*t - 0.09*t^2 + 0.001*t^3 + 50*ut
x <-  0.05 - 0.005*t + 0.09*t^2 - 0.001*t^3 + 40*ut
z <-  0.08 - 0.006*t + 0.08*t^2 - 0.001*t^3 + 30*ut

# The trend matrix
Trend <- cbind(1,poly(t,3,raw=TRUE))

# Estimating the model
stdlrm <- StDLM(y,cbind(x,z),lag=1,Trend=Trend,v=5,maxiter=2000)

# Generate arbitrary dates
dates <- seq(as.Date("2014/1/1"), as.Date("2016/1/1"), "weeks")

## Plotting the variable y, its estimated trend and the fitted value. 
lag <- 1
oldpar <- par(mfcol=c(3,1))
matplot(dates[(lag+1):length(y)],cbind(y[(lag+1):length(y)],stdlrm$fit,stdlrm$trend),
xlab="Months",type='l',lty=c(1,2,3),lwd=c(1,1,3),col=c("black","blue","black"),ylab=" ",xaxt="n")
axis.Date(1, at = seq(as.Date("2014/1/1"), as.Date("2016/1/1"), "months"),labels=TRUE)
legend("bottomleft",legend=c("data","trend","fitted values"),lty=c(1,2,3),lwd=c(1,1,3),
col=c("black","blue","black"),cex=.85)
hist(stdlrm$res,main="",xlab="") ## Histogram of y
matplot(dates[2:length(y)],cbind(stdlrm$cvar),xlab="Months",type='l',lty=2,lwd=1,
ylab="fitted variance",xaxt="n")
axis.Date(1, at = seq(as.Date("2014/1/1"), as.Date("2016/1/1"), "months"),labels=TRUE)
par(oldpar) 
}