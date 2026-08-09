\name{model2rjfun}
\alias{model2rjfun}
\alias{model2ssgrfun}
\alias{modelexpr}
\alias{SSmod2rjfun}
\title{model2rjfun}
\usage{
model2rjfun(modelformula, pvec, data = NULL, jacobian = TRUE, testresult = TRUE, ...)
SSmod2rjfun(modelformula, pvec, data = NULL, jacobian = TRUE, testresult = TRUE, ...)
model2ssgrfun(modelformula, pvec, data = NULL, gradient = TRUE, 
              testresult = TRUE, ...)
modelexpr(fun)
}
\arguments{
\item{modelformula}{A formula describing a nonlinear regression model.}

\item{pvec}{A vector of parameters.}

\item{data}{A dataframe, list or environment holding data used in the calculation.}

\item{jacobian}{Whether to compute the Jacobian matrix.}

\item{testresult}{Whether to test the function by evaluating it at \code{pvec}.}

\item{gradient}{Whether to compute the gradient vector.}

\item{fun}{A function produced by one of \code{model2rjfun} or \code{model2ssgrfun}.}

\item{\dots}{Dot arguments, that is, arguments that may be supplied by \code{name = value}
  to supply information needed to compute specific quantities in the model.}
}
  
\details{
If \code{pvec} does not have names, the parameters will have names
generated in the form \samp{p_<n>}, e.g. \code{p_1, p_2}.  Names that appear in
\code{pvec} will be taken to be parameters of the model.  

The \code{data} argument may be a dataframe, list or environment, or \code{NULL}.
If it is not an environment, one will be constructed using the components
of \code{data} with parent environment set to be
the environment of \code{modelformula}.  

\code{SSmod2rjfun} returns a function with header \code{function(prm)}, which
evaluates the residuals (and if \code{jacobian} is \code{TRUE} the
Jacobian matrix) of the selfStart model (the rhs is used) at \code{prm}.  
The residuals are defined to be the right hand side of \code{modelformula} 
minus the left hand side. Note that the selfStart model used in the model
formula must be available (i.e., loaded). If this function is called from 
\code{nlxb()} then the \code{control} element \code{japprox} must be 
set to value \code{SSJac}. 
}

\description{
These functions create functions to evaluate residuals or sums of squares at 
particular parameter locations.
}

\value{
\code{model2rjfun} returns a function with header \code{function(prm)}, which
evaluates the residuals (and if \code{jacobian} is \code{TRUE} the
Jacobian matrix) of the model at \code{prm}.  The residuals are defined to be
the right hand side of \code{modelformula} minus the left hand side.

\code{model2ssgrfun} returns a function with header \code{function(prm)}, which
evaluates the sum of squared residuals (and if \code{gradient} is \code{TRUE} the
gradient vector) of the model at \code{prm}. 

\code{modelexpr} returns the expression used to calculate the vector of 
residuals (and possibly the Jacobian) used in the previous functions.
}

\keyword{
   nonlinear
}

\examples{
# We do not appear to have an example for modelexpr. See nlsr-devdoc.Rmd for one.

y <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 38.558,
       50.156, 62.948, 75.995, 91.972)

tt <- seq_along(y)  # for testing

mydata <- data.frame(y = y, tt = tt)
f <- y ~ b1/(1 + b2 * exp(-1 * b3 * tt))
p <- c(b1 = 1, b2 = 1, b3 = 1)
rjfn <- model2rjfun(f, p, data = mydata)
rjfn(p)
rjfnnoj <- model2rjfun(f, p, data = mydata, jacobian=FALSE)
rjfnnoj(p)

myexp <- modelexpr(rjfn)
cat("myexp:"); print(myexp)

ssgrfn <- model2ssgrfun(f, p, data = mydata)
ssgrfn(p)

ssgrfnnoj <- model2ssgrfun(f, p, data = mydata, gradient=FALSE)
ssgrfnnoj(p)
}
\seealso{
\code{\link{nls}}
}
\author{
John Nash and Duncan Murdoch
}
