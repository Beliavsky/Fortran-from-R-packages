nlfb <-function(start, resfn, jacfn = NULL, trace = FALSE, 
		lower = -Inf, upper = Inf, weights=NULL,
		data=NULL, ctrlcopy=FALSE, control=list(), ...){
  
getweights <- function() { # local function
    if (is.null(weights0))
      rep(1, length(resraw))
    else if (is.function(weights0))
      weights0(pnum, resraw)
    else
      weights0
} 
  weights0 <- weights # save the original, possibly a function
  
  if (ctrlcopy) ctrl<-control else ctrl <- nlsr.control(control) 
  if (trace && ctrl$prtlvl>0) cat("psi=",ctrl$psi,"  phi=",ctrl$phi,"\n")
  # Function to display SS and point
  roff <- NA # initially undefined
  showprms<-function(SS, roff, pnum){
    pnames<-names(pnum)
#    npar<-length(pnum) # Will already be set before showprms called
    cat("lamda:",ctrl$lamda," SS=",SS,"(",roff,") at")
    for (i in 1:npar){
       cat(" ",pnames[i],"=",pnum[i])
    }
    cat(" f/j ",feval,"/",jeval)
    cat("\n")
  } # end showprms

  if (trace && ctrl$prtlvl>1) {
    if (is.null(weights)) {
       cat("no weights\n") 
    } else {
      cat("weights:")
      print(weights)
    }
  }
  # ensure params in vector
  if (missing(start)) stop("MUST provide starting parameters")
  pnames<-names(start)
  npar<-length(start)
  if (is.null(pnames)) {
    for (j in 1:npar) {  pnames[[j]]<- paste("p",j,sep='')}
  } 
  if (trace && ctrl$prtlvl>2) {cat("pnames:"); print(pnames)}
  pnum<-start # may simplify later
  names(pnum)<-pnames
  start<-as.numeric(start)
  names(start)<-pnames # needed?
  # bounds
  npar<-length(start) # number of parameters
  resraw <- resfn(pnum, ...)
  feval<-1 # number of function evaluations
  jeval<-0 # number of Jacobian evaluations
  if (length(lower)==1) lower<-rep(lower,npar)
  if (length(upper)==1) upper<-rep(upper,npar) 
  if (length(lower)!=npar) stop("Wrong length: lower")
  if (length(upper)!=npar) stop("Wrong length: upper")
  if (any(start < lower) || any(start > upper)) {
     if(trace) {  cat("start:"); print(start)
                     cat("\n");cat("lower:"); print(lower)
                     cat("\n"); cat("upper:"); print(upper) 
                  }
     stop("Infeasible start")
  }
  if (any(start<lower) || any(start>upper)) stop("Infeasible start")
  maskidx <- which(lower == upper) # MAY want to add tolerance?
  if (trace && ctrl$prtlvl>2) {cat("lower:"); print(lower); cat("upper:"); print(upper)}
  epstol<-(.Machine$double.eps)*ctrl$offset
  #   cat("ctrl$phi=",ctrl$phi,"\n")
  if (ctrl$phi == 0.0) phiroot<-0.0 else phiroot<-sqrt(ctrl$phi)
  if (ctrl$psi == 0.0) psiroot<-0.0 else psiroot<-sqrt(ctrl$psi)
  ctrl$stepredn <- ctrl$stepredn
  nbtlim <- ctrl$nbtlim 
  if (ctrl$stepredn <= 0) {
     nbtlim <- 1 # only 1 step if not backtracking
     if (trace) cat("No backtrack\n")
  } else {
     if (trace) cat("Bactrack with step reduction=",ctrl$stepredn,"\n")
  }
  bdmsk<-rep(1,npar) # set all params free for now
  if (length(maskidx)>0 && trace) {
       cat("The following parameters are masked:")
       print(pnames[maskidx]) # Change 140716
  }
  bdmsk[maskidx]<-0 # fixed parameters
  # Jacobian computations and approximations
  if (trace && ctrl$prtlvl>2) {cat("jacfn=",as.character(substitute(jacfn)),"\n") }
  appjac<-NULL # start NULL and test at end
  if (missing(jacfn) || is.null(jacfn)){
    if(is.null(ctrl$japprox)) stop("MUST PROVIDE jacobian function or specify approximation")
    if (trace && ctrl$prtlvl>0) cat("Using default jacobian approximation ",ctrl$japprox,"\n")
    if (trace) cat("ctrl$japprox=",ctrl$japprox,"\n")
    myjac<-match.fun(ctrl$japprox)
    appjac<-TRUE
  } else { 
## Need to see if jacfn has quotes and act appropriately
    if (is.character(jacfn)) { # have quoted jacobian function -- an approximation
      appjac <- TRUE # to inform us that we are using approximation
      myjac <- match.fun(jacfn)
      if (trace && ctrl$prtlvl>2) {cat("myjac="); print(myjac)}
    } else { # non-null, non-quoted jacfn
      if (trace && ctrl$prtlvl>2) {cat("We are using specified jacfn =",as.character(substitute(jacfn)),"\n")} 
      appjac<-FALSE # probably NOT needed
      myjac<-jacfn
    }
  }
  if (is.null(appjac)) stop("Failed to define jacfn -- appjac is null")
  if (trace && ctrl$prtlvl>1) {cat("Starting pnum="); print(pnum)}
  weights <- getweights()
  swts <- sqrt(weights)
  resbest<-resfn(pnum, ...)*swts
  ssbest<-as.numeric(crossprod(resbest))
  if (trace) { cat("00"); showprms(ssbest, roff, pnum) }
  ssminval <- ssbest*epstol^4
#  if (ctrl$watch) cat("ssminval =",ssminval,"\n")
  pbest<-pnum
  if (trace && ctrl$prtlvl>1) { 
    cat("Start:"); showprms(ssbest,roff,pnum);
    # if (ctrl$watch) tmp<-readline("Continue")
  }
  if (length(maskidx) == npar) {
    warning("All parameters are masked")
    result <- list(resid = resbest, jacobian = NA, feval = feval, 
            jeval = jeval, coefficients = pnum, ssquares = ssbest, 
            lower=lower, upper=upper, maskidx=maskidx, weights=weights0, 
            formula=NULL) 
    class(result)<-"nlsr"
    return(result)
  }
  ssquares<-.Machine$double.xmax # make it big
  newjac<-TRUE # set control for computing Jacobian
  eqcount<-0
  roffstop <- FALSE
  smallstop <- FALSE
  while ((! roffstop) && (eqcount < npar) && (feval <= ctrl$femax) 
           && (jeval <= ctrl$jemax) && (! smallstop) ) {
    ## if (trace) cat("Inside top main loop -- newjac=",newjac,"\n")
    if (newjac) {
      bdmsk<-rep(1,npar)
      bdmsk[which(pnum-lower<epstol*(abs(lower)+epstol))]<- -3 
      bdmsk[which(upper-pnum<epstol*(abs(upper)+epstol))]<- -1
      bdmsk[maskidx]<-0 # MUST be last -- or don't get masking
      if (ctrl$watch) { cat("bdmsk:"); print(bdmsk) }
      if (appjac) { # use approximation
       Jac<-myjac(pbest, resfn=resfn, bdmsk=bdmsk, 
                resbest=resbest, ndstep=ctrl$ndstep, ...)
      }
      else {Jac <- attr(myjac(pbest, ...),"gradient") }
      if (is.null(Jac)) stop("nlfb: Jac is not defined! (Is attr 'gradient' set?)")
      weights<-getweights()
      swts <- sqrt(weights)
      Jac <- Jac * swts
      resw <- resbest # already wtd
      ## NOTE: by insisting on using the "gradient" attribute, we can use same
      ## fn for gradient and residual
      jeval<-jeval+1 # count Jacobians
      if (any(is.na(Jac))) stop("NaN in Jacobian")
      gjty<-t(Jac)%*%resw # raw gradient
      if (trace && ctrl$prtlvl>2) { cat("gjty:"); print(as.vector(gjty)) }
      for (i in 1:npar){
        bmi<-bdmsk[i]
        if (bmi==0) {
          gjty[i]<-0 # masked
          Jac[,i]<-0
        }
        if (bmi<0) {
          if((2+bmi)*gjty[i] > 0) { # free parameter
            bdmsk[i]<-1
            if (ctrl$watch) cat("freeing parameter ",i," now at ",pnum[i],"\n")
          } 
          else {
            gjty[i]<-0 # active bound
            Jac[,i]<-0
            if (ctrl$watch) cat("active bound ",i," at ",pnum[i],"\n") 
          }
        } # bmi
      } # end for loop
      if (trace && ctrl$prtlvl>2) { cat("gjty-adj:"); print(as.vector(gjty)) }
      if (psiroot > 0.0) {
        if (npar == 1) dee <- diag(as.matrix(sqrt(diag(crossprod(Jac)))))
        else dee <- diag(sqrt(diag(crossprod(Jac))))  # to append to Jacobian
      }
    } # end newjac
    if (ctrl$jemax > 0) { # DO NOT COMPUTE CHANGE IF EVAL ONLY OPTION    
        lamroot<-sqrt(ctrl$lamda)
        JJ <- Jac
        rplus <- resw
        if (psiroot > 0.0) {
          JJ<-rbind(JJ,lamroot*psiroot*dee) # start building the matrix
          rplus <- c(rplus, rep(0, npar)) 
        }
        if (phiroot > 0.0) {
          JJ<-rbind(JJ,lamroot*phiroot*diag(npar)) # start building the matrix
          rplus <- c(rplus, rep(0, npar)) 
        }
        if (trace && ctrl$prtlvl>1) {
             cat("JJ\n"); print(JJ)
        }
        JQR<-qr(JJ)# Should we use try()?
        roff <- max(abs(as.numeric(crossprod(qr.Q(JQR), rplus))))/sqrt(ssbest+ctrl$scaleOffset)
        if (ctrl$watch) cat("roff =", roff,"  converged = ",(roff <= sqrt(epstol)),"\n")
        if (ctrl$rofftest && (roff <= sqrt(epstol))) roffstop <- TRUE
#       if (trace) cat("roffstop now ",roffstop,"\n")
        if (roffstop) break # added 220619
        #        tmp <- readline('cont')
        delta<-try(qr.coef(JQR,-rplus)) # Note appended rows of y)
        if (trace && ctrl$prtlvl>2) { cat("delta:"); print(delta) }
        if (inherits(delta, "try-error")) {
          if (ctrl$lamda<1000*.Machine$double.eps) ctrl$lamda<-1000*.Machine$double.eps
          ctrl$lamda<-ctrl$laminc*ctrl$lamda
          newjac<-FALSE # increasing lamda -- don't re-evaluate
          if (trace) cat(" Equation solve failure\n")
        } else { # solution OK
          gproj<-crossprod(delta,gjty)
          gangle <- gproj/sqrt(crossprod(gjty) * crossprod(delta))
          gangle <- 180 * acos(sign(gangle)*min(1, abs(gangle)))/pi
          if (ctrl$watch) cat("gradient projection = ",gproj," g-delta-angle=",gangle,"\n")
          if (is.na(gproj) || (gproj >= 0) ) { # uphill direction -- should NOT be possible
            if (ctrl$lamda<1000*.Machine$double.eps) ctrl$lamda <- 1000*.Machine$double.eps
            ctrl$lamda <- ctrl$laminc*ctrl$lamda
            newjac<-FALSE # increasing lamda -- don't re-evaluate
            if (trace) cat(" Uphill search direction\n")
            # May have trouble in bounded case!?
          } else { # downhill  
            delta[maskidx]<-0
            delta<-as.numeric(delta)
            if (ctrl$watch) {cat("delta:"); print(delta)}
            step<-rep(1,npar)
            for (i in 1:npar){
              bd<-bdmsk[i]
              da<-delta[i]
              #  if (ctrl$watch) cat(i," bdmsk=",bd,"  delta=",da,"\n")
              if (bd==0 || ((bd==-3) && (da<0)) ||((bd==-1)&& (da>0))) {
                delta[i]<-0
              } else {
                if (delta[i]>0) step[i]<-(upper[i]-pbest[i])/delta[i]
                if (delta[i]<0) step[i]<-(lower[i]-pbest[i])/delta[i] # positive
              }
            } # end loop over parameters/bounds
            stepsize<-min(1,step[which(delta != 0)]) # reduce stepsize if needed
            if (stepsize < .Machine$double.eps) {
               if (ctrl$lamda<1000*.Machine$double.eps) lamda<-1000*.Machine$double.eps
               ctrl$lamda<-ctrl$laminc*ctrl$lamda
               newjac<-FALSE # increasing lamda -- don't re-evaluate
               if (trace) cat(" Step too small or not possible; lamda now ",ctrl$lamda,"\n")
            }
            else { # stepsize not too small
              nbt=0 # number of backtracks (always do 1 evaluation)
              ssquares<-2*ssbest+1.0 # guarantee bigger
              while ((nbt < nbtlim) && (ssquares >= ssbest) ) { # backtrack loop 
                ## don't loop if ctrl$stepredn==0, and only until ssquares < ssbest
                pnum<-pbest+stepsize*delta # adjust step
                eqcount<-length(which((ctrl$offset+pbest)==(ctrl$offset+pnum))) 
                if (eqcount < npar) { # NOT all equal, so some change, computer resfn
                  feval<-feval+1 # count evaluations
                  nbt <- nbt+1
                  resid <- resfn(pnum, ...) * swts
                  ssquares<-sum(resid^2L) # get sum of squares
                  if (is.na(ssquares)) ssquares<-.Machine$double.xmax # Caution!
                  if (trace && ctrl$prtlvl>0)  cat(nbt," stepsize=",stepsize," lamda=",
                                   ctrl$lamda," ssquares=",ssquares,"\n")
                  if (ssquares >= ssbest) { # failed step, decrease stepsize
                     if ( (ctrl$stepredn <= 0.0) || (nbt >= ctrl$nbtlim) ) break # out of backtrack loop     
                     stepsize <- stepsize * ctrl$stepredn 
                  }
                } # end equality check
                else { # no change -- Do we need to fix in bounded case?
                  if (trace && ctrl$prtlvl>0) cat("No change in parameters(1)\n")
                  break # no change (eqcount == npar)
                }
              } # end while ssquares not reduced 
              # at this point we have either smaller ssquares (good), or failed step
              # failed step is either no change (converged) or we increase lambda
              if (ssquares < ssbest) { # smaller sumsquares, SUCCESS. Clean up and new iteration
                ctrl$lamda<-ctrl$lamdec*ctrl$lamda/ctrl$laminc  
                   # reduction via lamdec/laminc e.g., 4/10
                if (trace) { cat("<<"); showprms(ssquares, roff, pnum) }
                ssbest<-ssquares # save "best" so far
                if (ctrl$smallsstest) { smallstop<-(ssbest <= ssminval) } # capture small ss
                resbest<-resid
                pbest<-pnum
                newjac<-TRUE # indicate to compute new gradient
              } # end reduced sumsquares
              else { # Converged or increasing lambda
                if (eqcount < npar) {
                  if (ctrl$lamda<1000*.Machine$double.eps) ctrl$lamda<-1000*.Machine$double.eps
                  ctrl$lamda<-ctrl$laminc*ctrl$lamda
                  newjac<-FALSE # increasing lamda -- don't re-evaluate
                  if(trace && ctrl$prtlvl>2) showprms(ssquares, roff, pnum)
                }
                else {# equcount >= npar (will be ==). No progress.
                   if (trace) cat("No parameter change - stop!\n")
                }
              } # converged or increasing lambda
            } # end else stepsize not too small
          } # downhill
        } # solution OK
        if (ctrl$watch) tmp<-readline("Cycle")
     }
  } # end main while loop 
  pnum<-as.vector(pnum)
  names(pnum) <- pnames
  result <- list(resid = resbest, jacobian = Jac, feval = feval, 
            jeval = jeval, coefficients = pnum, ssquares = ssbest, lower=lower, 
            upper = upper, maskidx = maskidx, weights0 = weights0, 
            weights = weights, formula = NULL) # chg 190805
  class(result) <- "nlsr" ## Needed for print method
  result
} # end nlfb
