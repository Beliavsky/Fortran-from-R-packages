#' Compute a Relative Survival Curve from an additive relative survival model
#'
#' Computes the predicted relative survival function for an additive relative
#' survival model fitted with maximum likelihood.
#'
#' @param formula a \code{rsadd} object (Implemented only for models fitted
#' with the codemax.lik (default) option.)
#' @param newdata a data frame with the same variable names as those that
#' appear in the \code{rsadd} formula.  a predicted curve for each individual
#' in this data frame shall be calculated
#' @return a \code{survfit} object; see the help on \code{survfit.object} for
#' details.  The \code{survfit} methods are used for \code{print}, \code{plot},
#' \code{lines}, and \code{points}.
#' @seealso \code{survfit}, \code{survexp}
#' @references Package. Pohar M., Stare J. (2006) "Relative survival analysis
#' in R." Computer Methods and Programs in Biomedicine, \bold{81}: 272--278
#' @keywords survival
#' @examples
#'
#' data(slopop)
#' data(rdata)
#' #fit a relative survival model
#' fit <- rsadd(Surv(time,cens)~sex+age+year,rmap=list(age=age*365.241),
#' 	ratetable=slopop,data=rdata,int=c(0:10,15))
#'
#' #calculate the predicted curve for a male individual, aged 65, diagnosed in 1982
#' d <- rs.surv.rsadd(fit,newdata=data.frame(sex=1,age=65,year=as.Date("1982-01-01")))
#' #plot the curve (will result in a step function since the baseline is assumed piecewise constant)
#' plot(d,xscale=365.241)
#'
#' #calculate the predicted survival curves for each individual in the data set
#' d <- rs.surv.rsadd(fit,newdata=rdata)
#' #calculate the average over all predicted survival curves
#' p.surv <- apply(d$surv,1,mean)
#' #plot the relative survival curve
#' plot(d$time/365.241,p.surv,type="b",ylim=c(0,1),xlab="Time",ylab="Relative survival")
#'
rs.surv.rsadd <- function (formula, newdata)
{

  call <- match.call()
  Terms <- terms(formula$formula)  #to rabis, ce je model mal bl smotan - as.factor ali splines ali svasta
  Terms <- delete.response(Terms)
  newdata <- model.frame(Terms,newdata)
  n <- formula$n
  if(formula$method=="max.lik"){
    nvar <- length(formula$coef) - length(formula$int)+1
    formula$coef <- formula$coef[1:nvar]
  }
  nvar <- length(formula$coef)
  nx <- nrow(newdata)
  nt <- length(formula$times)
  temp <- list(n=formula$n,time=formula$times,call=call,type="right")
  Lambda0 <- formula$Lambda0
  Lambda0 <- matrix(Lambda0,ncol=nt,nrow=nx,byrow=TRUE)
  rate <- attr(Terms, "specials")$ratetable
  R <- as.matrix(newdata[, rate,drop=FALSE])
  rat <- attributes(formula$ratetable)$dimid
  mein <- attributes(newdata[,rate])$dimnames[[2]]
  x <- match(rat,mein)
  R <- R[,x,drop=FALSE]

  #newdata <- newdata[,1:nvar,drop=FALSE]
  labeli <- attr(attr(newdata,"terms"),"term.labels")
  colnami <- colnames(newdata)
  if(length(rate>0)){
    labeli <- labeli[-rate]
    colnami <- colnami[-rate]
  }
  newdata <- newdata[,match(colnami,labeli),drop=F]

  if(any(formula$mvalue)>0)newdata <- newdata - matrix(formula$mvalue,nrow=nx,byrow=TRUE)
  R <- data.frame(R)
  names(R) <- rat

  formula_covs <- formula(call("~", formula$formula[[length(formula$formula)]]))
  newdata2 <- model.matrix(formula_covs, newdata)
  newdata2 <- newdata2[, 2:ncol(newdata2), drop=FALSE]

  #ebx <- exp(data.matrix(newdata)%*%as.vector(formula$coef))
  ebx <- exp(data.matrix(newdata2)%*%as.vector(formula$coef))
  ebx <- matrix(ebx,ncol=nt,nrow=length(ebx))
  Lambdae <- Lambda0*ebx
  temp$surv <- t(exp(-Lambdae))
  temp$n.event <- rep(1,nt)
  temp$n.risk <- n+1 - cumsum(temp$n.event)
  temp$time <- formula$times
  class(temp) <- c("rs.surv.rsadd", "rs.surv","survfit")
  temp
}
