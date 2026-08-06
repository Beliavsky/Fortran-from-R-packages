#####################################################
## This file contains S3 methods for objects
## of class "tvgarchTest":
##
## coef.tvgarchTest()
## fitted.tvgarchTest()
## logLik.tvgarchTest()
## nobs.tvgarchTest()
## plot.tvgarchTest()
## predict.tvgarchTest()
## print.tvgarchTest()
## quantile.tvgarchTest()
## residuals.tvgarchTest()
## summary.tvgarchTest()
## toLatex.tvgarchTest()
## vcov.tvgarchTest() 
##
#####################################################


#####################################################
coef.tvgarchTest <- function (object, ...)
{
  return(coef.tvgarch(object = object$garch11))
}

#####################################################
fitted.tvgarchTest <- function (object, as.zoo = TRUE, ...)
{
  return(fitted.tvgarch(object = object$garch11, as.zoo = as.zoo))
}

#####################################################
logLik.tvgarchTest <- function(object, ...)
{
  return(logLik.tvgarch(object = object$garch11))
}

#####################################################
nobs.tvgarchTest <- function (object, ...)
{
  return(nobs.tvgarch(object = object$garch11))
}

#####################################################
plot.tvgarchTest <- function (x, ...) 
{
  plot.tvgarch(x = x$garch11)
}

#####################################################
predict.tvgarchTest <- function (object, n.ahead = 10, newxreg = NULL, 
                                 newindex = NULL, n.sim = 5000, as.zoo = TRUE, 
                                 verbose = FALSE, ...)
{
  return(predict.tvgarch(object = object$garch11, n.ahead = n.ahead, 
                         newxreg = newxreg, newindex = newindex, n.sim = n.sim, 
                         as.zoo = as.zoo, verbose = verbose))
}

#####################################################
print.tvgarchTest <- function (x, ...) 
{
  cat("\n")
  cat("Date:", x$date, "")
  cat("\n")
  cat("Testing GARCH(1,1), i.e., the model under H0, against 
  TV-GARCH(1,1), i.e., the model under H1: \n")
  cat("\n")
  cat("Estimation results for model under H0:")
  summary.tvgarch(object = x$garch11)
  cat("\n")
  cat("Transition variable in TV-GARCH(1,1):", colnames(x$xtv), "\n")
  cat("\n")
  cat("Results from the Non-Robust TR^2 Test: \n\n")
  print(x$mat)
  cat("\n")
  cat("Results from the Robust TR^2 Test: \n\n")
  print(x$matRob)
  cat("\n")
  print(x$order.g)
  invisible(x$order.g)
}

#####################################################
quantile.tvgarchTest <- function (x, probs = 0.025, names = TRUE, type = 7,
                                  as.zoo = TRUE, ...)
{
  return(quantile.tvgarch(x = x$garch11, probs = probs, names = names, 
                          type = type, as.zoo = as.zoo))
}

#####################################################
residuals.tvgarchTest <- function (object, as.zoo = TRUE, ...)
{
  return(residuals.tvgarch(object = object$garch11, as.zoo = as.zoo))
}

#####################################################
summary.tvgarchTest <- function (object, ...)
{
  cat("\n")
  cat("Date:", object$date, "")
  cat("\n")
  cat("Testing GARCH(1,1), i.e., the model under H0, against 
  TV-GARCH(1,1), i.e., the model under H1: \n")
  cat("\n")
  print(object$order.g)
  invisible(object$order.g)
}

#####################################################
toLatex.tvgarchTest <- function (object, digits = 4, ...)
{
  return(toLatex.tvgarch(object = object$garch11, digits = digits))
}

#####################################################
vcov.tvgarchTest <- function (object, ...)
{
  return(vcov.tvgarch(object = object$garch11))
}

