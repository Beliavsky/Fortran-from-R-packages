#' Compute a Predicited Survival Curve
#'
#' Computes a predicted survival curve based on the additive model estimated by
#' rsadd function.
#'
#' When predicting the survival curve, the ratetable values for future years
#' will be equal to those of the last given year. The same ratetables will be
#' used for fitting and predicting. To predict a relative survival curve, use
#' \code{rs.surv.rsadd}.
#'
#' @param formula a rsadd object
#' @param newdata a data frame with the same variable names as those that
#' appear in the rsadd formula. The curve(s) produced will be representative of
#' a cohort who's covariates correspond to the values in newdata.
#' @param se.fit a logical value indicating whether standard errors should be
#' computed. Default is \code{TRUE}.
#' @param conf.int the level for a two-sided confidence interval on the
#' survival curve(s). Default is 0.95.
#' @param individual a logical value indicating whether the data frame
#' represents different time epochs for only one individual (T), or whether
#' multiple rows indicate multiple individuals (F, the default). If the former
#' only one curve will be produced; if the latter there will be one curve per
#' row in newdata.
#' @param conf.type One of \code{none}, \code{plain}, \code{log} (the default),
#' or \code{log-log}. The first option causes confidence intervals not to be
#' generated. The second causes the standard intervals curve +- k *se(curve),
#' where k is determined from conf.int. The log option calculates intervals
#' based on the cumulative hazard or log(survival). The last option bases
#' intervals on the log hazard or log(-log(survival)).
#' @param ... Currently not implemented
#' @return a \code{survfit} object; see the help on \code{survfit.object} for
#' details.  The \code{survfit} methods are used for \code{print}, \code{plot},
#' \code{lines}, and \code{points}.
#' @seealso \code{survfit}, \code{survexp}, \code{\link{rs.surv}}
#' @references Package: Pohar M., Stare J. (2006) "Relative survival analysis
#' in R." Computer Methods and Programs in Biomedicine,\bold{81}: 272--278.
#'
#' Relative survival: Pohar, M., Stare, J. (2007) "Making relative survival
#' analysis relatively easy." Computers in biology and medicine, \bold{37}:
#' 1741--1749.
#' @keywords survival
#' @examples
#'
#' data(slopop)
#' data(rdata)
#' # work on a smaller dataset here to run the example faster
#' fit <- rsadd(Surv(time,cens)~sex,rmap=list(age=age*365.241),
#' 	ratetable=slopop,data=rdata[1:500,],method="EM")
#' survfit.rsadd(fit,newdata=rdata[1:2,])
#' survfit.rsadd(fit,newdata=data.frame(sex=1,age=60,year=17000))
#'
#'
survfit.rsadd <-
  function (formula, newdata, se.fit = TRUE, conf.int = 0.95, individual = FALSE,
            conf.type = c("log", "log-log", "plain", "none"),...)
  {

    call <- match.call()
    Terms <- terms(formula)  #to rabis, ce je model mal bl smotan - as.factor ali splines ali svasta
    Terms <- delete.response(Terms)
    popdata <- newdata
    newdata <- model.frame(Terms,newdata)
    resp <- list(y=formula$y,x=newdata)
    n <- formula$n
    nvar <- length(formula$coef)
    nx <- nrow(newdata)
    nt <- length(formula$times)
    temp <- list(n=formula$n,time=formula$times,call=call,type="right")
    Lambda0 <- formula$Lambda0
    Lambda0 <- matrix(Lambda0,ncol=nt,nrow=nrow(newdata),byrow=TRUE)

    rate <- attr(Terms, "specials")$ratetable

    #rat <- attributes(formula$ratetable)$dimid
    rat <- names(attributes(formula$ratetable)$dimnames)
    #mein <- attributes(newdata[,rate])$dimnames[[2]]
    mein <- names(popdata)
    x <- match(rat,mein)
    #R <- as.matrix(newdata[, rate, drop = FALSE])
    R <- popdata # <- as.matrix(popdata)
    R <- R[,x,drop=FALSE]
    # R <- data.frame(R)
    names(R) <- rat


    #newdata <- newdata[,1:(rate-1),drop=FALSE]
    labeli <- attr(attr(newdata,"terms"),"term.labels")
    colnami <- colnames(newdata)
    if(length(rate>0)){
      labeli <- labeli[-rate]
      colnami <- colnami[-rate]
    }
    newdata <- newdata[,match(colnami,labeli),drop=F]
    if(any(formula$mvalue)>0)newdata <- newdata - matrix(formula$mvalue,nrow=nrow(newdata),byrow=TRUE)
    nx <- ncol(newdata)
    #getl <- function(times,data=R,ratetable=formula$ratetable){
    #	-log(srvxp.fit(data,times,ratetable))
    #}
    #Lambdap <- sapply(formula$times, getl)
    # Lambdap <- NULL
    # for(it in 1:nt){
    # 	Lambdap <- cbind(Lambdap,-log(srvxp.fit(R,formula$times[it],formula$ratetable)))
    # }

    Lambdap <- NULL
    for(it in 1:nrow(newdata)){
      Lambdap <- rbind(Lambdap,-log(survexp(~1,data=R[it,,drop=FALSE],times=formula$times,ratetable=formula$ratetable)$surv))
    }
    # for(ie in 1:ncol(newdata)){
    #   newdata[,ie] <- as.numeric(newdata[,ie])
    # }
    formula_covs <- formula(call("~", formula$formula[[length(formula$formula)]]))
    newdata2 <- model.matrix(formula_covs, newdata)
    newdata2 <- newdata2[, 2:ncol(newdata2), drop=FALSE]

    # ebx <- exp(as.matrix(data.frame(coef=formula$coef)^T %*%newdata))
    ebx <- exp(newdata2 %*% matrix(formula$coef))
    ebx <- matrix(ebx,ncol=nt,nrow=length(ebx))
    Lambda <- Lambdap + Lambda0*ebx
    temp$surv <- t(exp(-Lambda))
    temp$n.event <- rep(1,nt)
    temp$n.risk <- n+1 - cumsum(temp$n.event)
    class(temp) <- c("rs.surv.rsadd", "rs.surv","survfit")
    temp
  }
