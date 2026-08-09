nlsr.package <- function(){
   cat("nlsr package: a set of tools for solving nonlinear least squares problems\n")
   return(1)
}

nlsr.control <- function(control) {
  # NEED TO PUT IN CHECKS of validity
  # Defaults
  defctrl=list(femax=10000,
              japprox=NULL,
              jemax=5000,
              lamda=0.0001,
	      laminc=10, 
	      lamdec=4, 
          nbtlim=6,
	      ndstep=1e-7,
	      offset=100, 
	      phi=1.0, 
	      psi=0.0,
	      rofftest = TRUE, 
	      smallsstest = TRUE, 
	      stepredn=0.0,
	      watch=FALSE,
              prtlvl=0,
              scaleOffset=1.0)
  if (missing(control) || is.null(control)){control <- defctrl}
  else {
      namctl <- names(control)
      namdef <- names(defctrl)
      if (! all(namctl %in% namdef) ) stop("unrecognized control names")
      for (onename in namctl) {defctrl[onename]<-control[onename]}
      control<-defctrl
  }
  control # return the control list
}

coef.nlsr <- function(object, ...) {
       out <- object$coefficients
       attr(out,"pkgname")<-"nlsr"
       out # JN 170109
}

summary.nlsr <- function(object, ...) {
  smalltol <- .Machine$double.eps * 1000
  options(digits = 5) # 7 is default
  resname <- deparse(substitute(object))
  JJ <- object$jacobian
  res <- object$resid # weighted
  coeff <- object$coefficients
  pnames<-names(coeff)
#  cat("pnames:"); print(pnames)
  npar <- length(coeff)
  w <- object$weights
  nobs <- if (!is.null(w)) sum(w > 0) else length(res)
  lo <- object$lower
  if (is.null(lo)) lo <- rep( -Inf, npar)
  up <- object$upper
  if (is.null(up)) up <- rep( Inf, npar)
  mi <- object$maskidx
  mt <- rep(" ",npar) # start with all "unmasked"
  mt[mi] <- "M" # Put in the masks
  bdmsk <- rep(1, npar) # bounds and masks indicator. Should it be 1L?
  bdmsk[mi] <- 0 # masked
  ct <- rep(" ",npar) # start with all "free"
  for (i in seq_along(coeff)){
    if (lo[[i]] - coeff[[i]] > 0) {
      ct[[i]] <- "-" # lower bound violation
      if (bdmsk[[i]] == 1) bdmsk[[i]] <- -3
    } else { 
      if (coeff[[i]] - lo[[i]] < smalltol*(abs(coeff[[i]])+smalltol) ) {
        ct[[i]] <- "L" # "at" lower bound
        if (bdmsk[[i]] != 0) bdmsk[[i]] <- -3 # leave mask indication intact
      }
    }
    if (coeff[[i]] - up[[i]] > 0) {
      ct[[i]] <- "+" # upper bound violation
      if (bdmsk[[i]] == 1) bdmsk[[i]] <- -1
    } else { 
      if (up[[i]] - coeff[[i]] < smalltol*(abs(coeff[[i]])+smalltol) ) {
        ct[[i]] <- "U" # "at" upper bound
        if (bdmsk[[i]] != 0) bdmsk[[i]] <- -1 # leave mask indication intact
      }
    }
  }
  notmask<-which(bdmsk != 0)
  ss <- object$ssquares
  rdf <- nobs - length(notmask) # New way 220801
  if (rdf <= 0) {
    if (rdf < 0) { stop(paste("Inadmissible degrees of freedom =",rdf,sep='')) }
    else { resvar <- Inf }
  } else {
    resvar <- ss/rdf
  }
  nanc<-rep(NA,npar)
  Sdout <- SEs <- nanc # Set so it is defined and right length. Would Inf be preferred?  
  XtXinv <- matrix(NA, nrow=npar, ncol=npar) # default. 
  gr <- nanc # default to NA
  if (length(notmask) < 1) {
    warning("No unmasked parameters")
    XTXinv <- NULL
    param <- cbind(coeff, nanc, nanc, nanc)
  }
  else {
    Jwork<-JJ[, notmask]
    if (any(is.na(Jwork))){ stop("Bad working Jacobian") }
    else {
      dec <- svd(Jwork)
      U <- dec$u
      V <- dec$v
      Sd <- dec$d
    }
#    if (min(Sd) <= smalltol * max(Sd)) { # singular
#      XtXinv <- matrix(NA, nrow=npar, ncol=npar)
#    } else {  ## Above not needed because we have set default
    if (min(Sd) > smalltol * max(Sd)) { # singular
        Sinv <- 1/Sd
      if (length(notmask) > 1)  { # 220809 use 
        VS <- crossprod(t(V), diag(Sinv))
      } else {
        VS <- V/Sinv
      }
      XtXinv <- crossprod(t(VS))
      SEs[notmask] <- sqrt(diag(XtXinv) * resvar) # SEs of masked params still NA
    }
    Sdout[notmask] <- Sd
    gr <- crossprod(JJ, res)
    tstat<-rep(NA, npar)
    tstat[notmask] <- coeff[notmask]/SEs[notmask]
    tval <- coeff/SEs
    if (rdf > 0) {
       pval <- pt(abs(tval), rdf, lower.tail = FALSE)
    } else { pval <- NaN }
    param <- cbind(coeff, SEs, tval, 2 * pval)
  }
  # Note: We don't return formula because we may be doing nlfb summary 
  #   i.e., resfn and jacfn approach.
#  dimnames(param) <-
#    list(pnames, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  colnames(param) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  rownames(param) <- pnames
  ans <- list(residuals = res, sigma = sqrt(resvar),  
              df = c(npar, rdf), cov.unscaled = XtXinv,
              param = param, resname=resname, ssquares=ss, nobs=nobs, 
              ct=ct, mt=mt, Sd=Sdout, gr=gr, jeval=object$jeval,feval=object$feval)
  class(ans)<-"nlsr"
  ans
} # end summary()

print.nlsr <- function(x, ...) {
  if (inherits(x,"try-error") || is.null(x$coefficients)) {
    cat("Object has try-error or missing parameters\n")
    return(invisible(x))
  }
  xx<-summary(x) # calls summary to get information
  with(xx, { 
    pname<-rownames(param) # param is augmented coefficients with SEs and tstats
    npar <- dim(param)[1] # previously length(coeff) 
    cat("residual sumsquares = ",ssquares," on ",nobs,"observations\n")
    cat("    after ",jeval,"   Jacobian and ",feval,"function evaluations\n")
    cat("  name     ","      coeff    ","     SE   ","   tstat  ",
        "   pval  ","   gradient  "," JSingval  ","\n")
    SEs <- param[,2]
    tstat <- param[,3]
    pval <- param[,4]
#    for (i in seq_along(param[,1])){
    for (i in 1:dim(param)[1]){
      tmpname<-pname[i]
      if (is.null(tmpname)) {tmpname <- paste("p_",i,sep='')}
      cat(format(tmpname, width=10)," ")
      cat(format(param[[i]], digits=6, width=12))
      cat(ct[[i]],mt[[i]]," ")
      cat(format(SEs[[i]], digits=4, width=9)," ")
      cat(format(tstat[[i]], digits=4, width=9)," ")
      cat(format(pval[[i]], digits=4, width=9)," ")
      cat(format(gr[[i]], digits=4, width=10)," ")
      cat(format(Sd[[i]], digits=4, width=10)," ")
      cat("\n")
    }
  }) # remember to close with()
  # return(NULL) 
  invisible(x)
} # end print.nlsr()

pshort <- function(x) {
  xname<-deparse(substitute(x))
  cat(xname," -- ss=",x$ssquares,":")
  prm<-x$coefficients
  pnam<-names(prm)
  for (i in 1:length(prm)){ cat(" ",pnam[i],"=",prm[i]) }
  cat(";",x$feval,"res/",x$jeval,"jac\n")
  invisible(x)
}


resid.nlsr <- function(object, ...){ # weighted
  resids <- object$resid
  resids # so the function prints
}

residuals.nlsr <- function(object, ...){
  resids <- object$resid
  resids # so the function prints
}


predict.nlsr <- function(object=NULL, newdata=list(), ...) { 
#  This ONLY works if we have used nlxb. Do we want to add class 'nlxb'?
    if( is.null(object) ) stop("predict.nlsr REQUIRES an nlsr solution object")
    form <- object$formula
    if (is.null(form)) stop("nlsr.predict works only if formula is defined")
# Give more output?
#  Note: we assume a formula of style y~something, and use the something
#  In some ways need to check this more carefully
    env4coefs <- list2env(as.list(object$coefficients))
    preds <- eval(form[[3]], as.list(newdata), env4coefs)
    class(preds)<- "predict.nlsr"
    attr(preds,"pkgname") <- "nlsr"
    preds
}

fitted.nlsr <- function(object=NULL, data=parent.frame(), ...) { 
  #  This ONLY works if we have used nlxb. Do we want to add class 'nlxb'?
  # Is parent.frame() right?  What about dotargs?
  if( is.null(object) ) stop("predict.nlsr REQUIRES an nlsr solution object")
  form <- object$formula
  if (is.null(form)) stop("fitted.nlsr works only if formula is defined")
  data <- eval(object$data)
  # Give more output
  #  we assume a formula of style y~something, and use the something
  env4coefs <- list2env(as.list(object$coefficients))
  fits <- eval(form[[3]], as.list(data), env4coefs)
  class(fits)<- "predict.nlsr" # Do we want all this?
#  attr(fits,"pkgname") <- "nlsr" # Do we want all this?
  fits
}

rawres <- function(object=NULL, data=parent.frame(), ...) { 
  #  This ONLY works if we have used nlxb. Do we want to add class 'nlxb'?
  # Is parent.frame() right?  What about dotargs?
  if( is.null(object) ) stop("predict.nlsr REQUIRES an nlsr solution object")
  form <- object$formula
  data <- eval(object$data) # Nov 17, 2022
  if (is.null(form)) stop("fitted.nlsr works only if formula is defined")
  env4coefs <- list2env(as.list(object$coefficients))
  rawres <- eval(form[[3]]-form[[1]], as.list(data), env4coefs)
  class(rawres)<- "nlxb" # Do we want all this?
  rawres
}

pctrl <- function(control) { 
   nprow <- 4 # default to 4 per row. May change later.
   nc <- length(control) # note some items may be multiple
   namctl<-names(control)
   for (ii in 1:nc) {
      item <- control[[ii]]
      if (length(item) > 1) item <- "vec/list"
      cat("  ",namctl[[ii]],"=",item)
      if ((ii %% nprow) == 0) cat("\n")
   }
   if ((ii %% nprow) != 0) cat("\n")
   invisible(1)
}

pnlslm <- function(x){ # short output for nls.lm result
  xname <- deparse(substitute(x))
  ss<-deviance(x)
  prm<-coef(x)
  pnam <- names(prm)   
  cat(xname," -- ss=",ss,":")
  for (i in 1:length(prm)){cat(" ",pnam[i],"=",prm[i])}
  cat(";",x$niter," itns\n") # only diff from pnlsLM
  invisible(1)
}

pnls <- function(x){ # short output for nls() result
  xname <- deparse(substitute(x))
  ss<-deviance(x)
  prm<-coef(x)
  pnam <- names(prm)   
  cat(xname," -- ss=",ss,":")
  for (i in 1:length(prm)){cat(" ",pnam[i],"=",prm[i])}
  cat(";",x$convInfo$finIter," itns\n")
  invisible(1)
}

# pnlsLM <- function(x){ # short output for nlsLM result
#   xname <- deparse(substitute(x))
#   ss<-deviance(x)
#   prm<-coef(x)
#   pnam <- names(prm)   
#   cat(xname," -- ss=",ss,":")
#   for (i in 1:length(prm)){cat(" ",pnam[i],"=",prm[i])}
#   cat(";",x$convInfo$finIter," itns\n")
#   invisible(1)
# }
# 


nvec <- function(vec){ # tidy display for named vector
   vnam<-deparse(substitute(vec))
   cat(vnam,":")
   pnam <- names(vec)
   for (i in 1:length(pnam)){ cat(pnam[i],"=",as.numeric(vec[i])," ")}
   cat("\n")
   invisible(1)
}

prt <- function(x) { # to provide naming for generic print 221123
  xname<-deparse(substitute(x))
  cat("nlsr object:",xname,"\n")
  print(x)
  invisible(x)
}
