jand <- function(pars, resfn=NULL, bdmsk=NULL, resbest=NULL, ndstep=1e-7, ...){
   # ndstep not currently used
   if (is.null(resfn)) stop("You must supply a residual function")
   # numDeriv jacobian approximation
   jacmat<-numDeriv::jacobian(resfn, pars,...) # don't specify other numDeriv parameters here yet
   attr(jacmat,"gradient") <- jacmat # to satisfy needs of nlxb(), nlfb() working
} # end jand

jafwd <- function(pars, resfn=NULL, bdmsk=NULL, resbest=NULL, ndstep=1e-7, ...){
  # Forward difference jacobian approximation
   if (is.null(resfn)) stop("You must supply a residual function")
   npar<-length(pars)
   nres<-length(resbest)
   jacmat<-matrix(0, ncol=npar, nrow=nres)
   tpars<-pars # need to set because values used for each dimension
   for (j in 1:npar){
      if (is.null(bdmsk) || (bdmsk[j]!=0)) {
         step<-ndstep*(abs(pars[j])+ndstep)
         tpars[j]<-pars[j]+step
         jacmat[,j]<-(resfn(tpars,...)-resbest)/step
         tpars[j]<-pars[j]
      }
   }
   attr(jacmat,"gradient") <- jacmat # to satisfy needs of nlxb(), nlfb() working
   jacmat
} # end jafwd

jaback <- function(pars, resfn=NULL, bdmsk=NULL, resbest=NULL, ndstep=1e-7, ...){
   # Backward difference jacobian approximation
   if (is.null(resfn)) stop("You must supply a residual function")
   npar<-length(pars)
   nres<-length(resbest)
   jacmat<-matrix(0, ncol=npar, nrow=nres)
   tpars<-pars
   for (j in 1:npar){
      if (is.null(bdmsk) || (bdmsk[j]!=0)) {
         step<-ndstep*(abs(pars[j])+ndstep)
         tpars[j]<-tpars[j]-step
         jacmat[,j]<-(resbest-resfn(tpars,...))/step
         tpars[j]<-pars[j]
      }
   }
   attr(jacmat,"gradient") <- jacmat # to satisfy needs of nlxb(), nlfb() working
   jacmat
} # end jaback

jacentral <- function(pars, resfn=NULL, bdmsk=NULL, resbest=NULL, ndstep=1e-7, ...){
   # Central difference jacobian approximation
   if (is.null(resfn)) stop("You must supply a residual function")
   npar<-length(pars)
   nres<-length(resbest)
   jacmat<-matrix(0, ncol=npar, nrow=length(resbest))
   tparf<-pars
   tparm<-pars
   for (j in 1:npar){
      if (is.null(bdmsk) || (bdmsk[j]!=0)) {
         step<-ndstep*(abs(pars[j])+ndstep)
         tparf[j]<-tparf[j]+step
         tparm[j]<-tparm[j]-step
         jacmat[,j]<-0.5*(resfn(tparf,...)-resfn(tparm,...))/step
         tparf[j]<-pars[j]
         tparm[j]<-pars[j]
      }
   }
   attr(jacmat,"gradient") <- jacmat # to satisfy needs of nlxb(), nlfb() working
   jacmat
} # end jacentral


resgr <- function(prm, resfn, jacfn, ...) {
    # computes the gradient 2 * J' %*% res for residuals (res)
    # and jacobian (Jac) defined by resfn and jacfn at
    # parameters prm and extra variables defined in the
    # dot-arguments J C Nash 2012-4-26
    res <- resfn(prm, ...)  # Should we use try()?
    if (missing(jacfn) || is.null(jacfn)){
       stop("MUST PROVIDE jacobian function or specify approximation")
    } else {  # Need to see if jacfn has quotes
      if (is.character(jacfn)) { # have quoted jacobian function -- an approximation
        myjac <- match.fun(jacfn)
        resval <- resfn(prm)
        cat("myjac:"); print(myjac)
        Jac <- myjac(prm, resfn, resbest=resval, ...)
      } else { # non-null, non-quoted jacfn
        Jac<-jacfn(prm, ...)
        cat("jacfn:"); print(jacfn)
      }
    }
    grj <- 2 * as.numeric(crossprod(Jac, res))
    attr(res, "Jacobian") <- Jac 
    attr(res, "gradient") <- grj
    res
}


resss <- function(prm, resfn, ...) {
    # computes sumsquares function from residuals defined by
    # resfn at parameters prm and extra variables defined in
    # the dot-arguments J C Nash 2012-4-26
    resids <- resfn(prm, ...)  #  try() ?
    ss <- as.numeric(crossprod(resids))
}

