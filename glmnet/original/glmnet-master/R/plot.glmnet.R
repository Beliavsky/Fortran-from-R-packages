#' plot coefficients from a "glmnet" object
#'
#' Produces a coefficient profile plot of the coefficient paths for a fitted
#' \code{"glmnet"} object.
#'
#' A coefficient profile plot is produced. If \code{x} is a multinomial model,
#' a coefficient plot is produced for each class.
#'
#' @aliases plot.glmnet plot.multnet plot.mrelnet plot.relaxed
#' @param x fitted \code{"glmnet"} model
#' @param xvar What is on the X-axis. \code{"lambda"} plots against the log-lambda sequence,
#' \code{"norm"} against the L1-norm of the coefficients, and
#' \code{"dev"} against the percent deviance explained. Warning: "norm" is the L1 norm of the coefficients on the glmnet object. There are many reasons why this might not be appropriate, such as automatic standardization, penalty factors, and values of `alpha` less than 1, which can lead to unusual looking plots.
#' @param label If \code{TRUE}, label the curves with variable sequence
#' numbers.
#' @param sign.lambda If \code{xvar="lambda"} and \code{sign.lambda=1} then we plot against \code{log(lambda)}; if  \code{sign.lambda=-1} (default) we plot against \code{-log(lambda)}.
#' @param \dots Other graphical parameters to plot
#' @author Jerome Friedman, Trevor Hastie and Rob Tibshirani\cr Maintainer:
#' Trevor Hastie <hastie@@stanford.edu>
#' @seealso \code{glmnet}, and \code{print}, \code{predict} and \code{coef}
#' methods.
#' @references Friedman, J., Hastie, T. and Tibshirani, R. (2008)
#' \emph{Regularization Paths for Generalized Linear Models via Coordinate
#' Descent}
#' @keywords models regression
#'
#' @examples
#' x=matrix(rnorm(100*20),100,20)
#' y=rnorm(100)
#' g2=sample(1:2,100,replace=TRUE)
#' g4=sample(1:4,100,replace=TRUE)
#' fit1=glmnet(x,y)
#' plot(fit1)
#' plot(fit1,xvar="lambda",label=TRUE)
#' fit3=glmnet(x,g4,family="multinomial")
#' plot(fit3,pch=19)
#' @method plot glmnet
#' @export
plot.glmnet=function(x, xvar=c("lambda","norm","dev"), label=FALSE, sign.lambda=-1,...){
  xvar=match.arg(xvar)
  plotCoef(x$beta,lambda=x$lambda,df=x$df,dev=x$dev.ratio,label=label,xvar=xvar,
           sign.lambda=sign.lambda,...)
}
