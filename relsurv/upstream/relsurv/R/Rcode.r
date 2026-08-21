rsfitterem<-function(data,b,maxiter,ratetable,tol,bwin,p,cause,Nie){
  # cause: = 2 (unknown), 0 in 1 known. Lahko preko argumenta cause v rsadd dolocis, ce kdo ima znan cause of death (ne rabijo vsi) .
  # Nie: to je lambda_0 (ti), ki se oceni v M koraku v EM algoritmu

  pr.time<-proc.time()[3]
  if (maxiter<1) stop("There must be at least one iteration run")
  n<-nrow(data)
  m <- p
  dtimes <- which(data$stat==1)			#the positions of event times in data$Y
  td <- data$Y[dtimes]				#event times
  ntd <- length(td)				#number of event times
  utimes <- which(c(1,diff(td))!=0)		#the positions of unique event times among td
  utd <- td[utimes]				#unique event times
  nutd <- length(utd)				#number of unique event times
  udtimes <- dtimes[utimes]			#the positions of unique event times among data$Y
  razteg <- function(x){
    # x is a 0/1 vector, the output is a vector of length sum(x), with the corresponding rep numbers
    n <- length(x)
    repu <- rep(1,n)
    repu[x==1] <- 0
    repu <- rev(cumsum(rev(repu)))
    repu <- repu[x==1]
    repu <- -diff(c(repu,0))+1
    if(sum(repu)!=n)repu <- c(n-sum(repu),repu)  #ce je prvi cas censoring, bo treba se kej narest??
    repu
  }
  rutd <- rep(0,ntd)
  rutd[utimes] <- 1
  # the number of times every unique event time occurs:
  rutd <- razteg(rutd)				#from  unique event times to event times
  # the number of censorings that occur after some event:
  rtd <- razteg(data$stat)			#from event times to data$Y

  a <- data$a[data$stat==1]

  if(bwin[1]!=0){
    #the vector of change points for the smoothing bandwidth
    nt4 <- c(1,ceiling(c(nutd*.25,nutd/2,nutd*.75,nutd)))
    if(missing(bwin))bwin <- rep(1,4)
    else bwin <- rep(bwin,4)
    for(it in 1:4){
      bwin[it] <- bwin[it]*max(diff(utd[nt4[it]:nt4[it+1]]))
    }
    while(utd[nt4[2]]<bwin[1]){		# ce je bwin velik, skrajsamo nt4
      nt4 <- nt4[-2]
      if(length(nt4)==1)break
    }
    #the smoothing matrix
    krn <- kernerleftch(utd,bwin,nt4)
  }



  #forming the new dataset
  if(p>0){
    whtemp <- data$stat==1&cause==2
    dataded <- data[data$stat==1&cause==2,]				#events with unknown cause
    datacens <- data[data$stat==0|cause<2,]				#censorings or known cause

    datacens$cause <- cause[data$stat==0|cause<2]*data$stat[data$stat==0|cause<2]

    databig <- lapply(dataded, rep, 2)
    databig <- do.call("data.frame", databig)
    databig$cause <- rep(2,nrow(databig))
    nded <- nrow(databig)
    databig$cens <- c(rep(1,nded/2),rep(0,nded/2))

    datacens$cens <- rep(0,nrow(datacens))
    datacens$cens[datacens$cause<2] <- datacens$cause[datacens$cause<2]
    names(datacens) <- names(databig)

    databig <- rbind(databig,datacens)

    cause <- cause[data$stat==1]

    #NEW IN 2.05 (next 4 lines)
    fk <- (attributes(ratetable)$factor != 1)
    nfk <- length(fk)
    varstart <- 3+nfk+1		#first column of covariates
    varstop <- 3+nfk+m		#last column of covariates
    #model matrix for relative survival
    xmat <- as.matrix(data[,varstart:varstop])		#NEW IN 2.05

    #ebx at initial values of b
    ebx <- as.vector(exp(xmat%*%b)) # exp(linear.predictor)

    #model matrix for coxph
    modmat <- as.matrix(databig[,varstart:varstop])		#NEW IN 2.05
    varnames <- names(data)[varstart:varstop]		#NEW IN 2.05
  }
  else{
    cause <- cause[data$stat==1]
    ebx <- rep(1,n) # exp(linear.predictor)
  }


  #for time-dependent data:
  starter <- sort(data$start)
  starter1<-c(starter[1],starter[-length(starter)])

  #the values of interest in the cumsums of the obsolete values (there is at least one value - the 1st)
  index <- c(TRUE,(starter!=starter1)[-1])

  starter <- starter[index]
  #the number of repetitions in each cumsum difference - needed for s0 calculation
  val1 <- apply(matrix(starter,ncol=1),1,function(x,Y)sum(x>=Y),data$Y)
  val1 <- c(val1[1],diff(val1),length(data$Y)-val1[length(val1)])

  eb <- ebx[data$stat==1]
  # s0 je sum_{at risk set} ebx
  s0 <- cumsum((ebx)[n:1])[n:1]

  ebx.st <- ebx[order(data$start)]
  s0.st <- ((cumsum(ebx.st[n:1]))[n:1])[index]
  s0.st <- rep(c(s0.st,0),val1)
  s0 <- s0 - s0.st


  #s0 only at times utd
  s0 <- s0[udtimes]


  #find the corresponding value of Y for each start!=0 - needed for likelihood calculation
  start <- data$start
  # if(any(start!=0)){
  #   wstart <- rep(NA,n)
  #   ustart <- unique(start[start!=0])
  #   for(its in ustart){
  #     wstart[start==its] <- min(which(data$Y==its))
  #   }
  # }


  #tale del je zelo sumljiv - kako se racuna likelihood za ties???
  difft <- c(data$Y[data$stat==1][1],diff(td))
  difft <- difftu <- difft[difft!=0]
  difft <- rep(difft,rutd)
  a0 <- a*difft



  if(sum(Nie==.5)!=0)maxit0 <- maxiter
  else maxit0<- maxiter - 3
  for(i in 1:maxit0){

    #Nie is of length ntd, should be nutd, with the values at times being the sum
    nietemp <- rep(1:nutd,rutd)
    Nies <- as.vector(by(Nie,nietemp,sum))  #shorter Nie - only at times utd

    lam0u <- lam0 <- Nies/s0
    #the smooting of lam0
    if(bwin[1]!=0){
      lam0s <- krn%*%lam0
    } else{
      lam0s <- lam0/difftu
    }

    #extended to all event times
    lam0s <- rep(lam0s,rutd)


    #compute Nie, only for those with unknown hazard
    Nie[cause==2] <- as.vector(lam0s*eb/(a+lam0s*eb))[cause==2]

  }

  if(maxit0!=maxiter & i==maxit0) i <- maxiter
  #likelihood calculation - manjka ti se likelihood za nicelni model!!!
  #the cumulative hazard
  Lam0  <- cumsum(lam0)
  #extended to all event times
  Lam0 <- rep(Lam0,rutd)
  if(data$stat[1]==0) Lam0 <- c(0,Lam0)
  #extended to all exit times
  Lam0 <- rep(Lam0,rtd)
  # save original object:
  Lam00 <- Lam0

  #for time dependent covariates and left-truncated individuals: replace by the difference
  if(any(start!=0)){
    # Calculate hazards at non-event times:
    timehaz <- data.frame(time=sort(data$Y), Lam0_2=Lam0)
    timehaz_tmp <- data.frame(time=unique(data$start), Lam0_2=NA)
    timehaz <- rbind(timehaz, timehaz_tmp)
    timehaz <- timehaz[order(timehaz$time),]
    timehaz$Lam0_2 <- mstateNAfix(timehaz$Lam0_2, 0)
    timehaz <- timehaz[!duplicated(timehaz$time),]

    # Prepare object so that you can calculate Lam0_event_time - Lam0_entry_time
    data_lt <- cbind(data, Lam0, id_0=1:nrow(data))
    data_lt <- merge(data_lt, timehaz,
                     by.x='start', by.y='time', all.x = TRUE)
    data_lt <- data_lt[order(data_lt$id_0),]

    # Check:
    # if(any(data_lt$Lam0_2[start!=0] != Lam0[wstart[start!=0]])){
    #   browser()
    # }

    # Edit Lam0:
    Lam0[start!=0] <- data_lt$Lam0[start!=0] - data_lt$Lam0_2[start!=0]

    # Old calculation:
    # Lam0[start!=0] <- Lam0[start!=0] - Lam0[wstart[start!=0]]
  }

  lam0 <- rep(lam0,rutd)

  likely0 <- sum(log(a0 + lam0*eb)) - sum(data$ds + Lam0*ebx)
  likely <- likely0
  tempind <- Nie<=0|Nie>=1
  if(any(tempind)){
    if(any(Nie<=0))Nie[Nie<=0] <- tol
    if(any(Nie>=1))Nie[Nie>=1] <- 1-tol
  }

  if(p>0)databig$wei <- c(Nie[cause==2],1-Nie[cause==2],rep(1,nrow(datacens)))


  if(maxiter>=1&p!=0){
    for(i in 1:maxiter){

      if(p>0){
        b00<-b
        if(i==1)fit <- coxph(Surv(start,Y,cens)~modmat,data=databig,weights=databig$wei,init=b00,x=TRUE,iter.max=maxiter)
        else    fit <- coxph(Surv(start,Y,cens)~modmat,data=databig,weights=databig$wei,x=TRUE,iter.max=maxiter)

        if(any(is.na(fit$coeff))) stop("X matrix deemed to be singular, variable ",which(is.na(fit$coeff)))

        b <- fit$coeff

        ebx <- as.vector(exp(xmat%*%b))
      }
      else ebx <- rep(1,n)

      eb <- ebx[data$stat==1]

      # s0 je sum_{at risk set} ebx
      s0 <- cumsum((ebx)[n:1])[n:1]

      ebx.st <- ebx[order(data$start)]
      s0.st <- ((cumsum(ebx.st[n:1]))[n:1])[index]
      s0.st <- rep(c(s0.st,0),val1)
      s0 <- s0 - s0.st

      #Nie is of length ntd, should be nutd, with the values at times being the sum
      nietemp <- rep(1:nutd,rutd)
      Nies <- as.vector(by(Nie,nietemp,sum))  #shorter Nie - only at times utd
      #s0 only at times utd
      s0 <- s0[udtimes]

      lam0u <- lam0 <- Nies/s0


      #the cumulative hazard
      Lam0  <- cumsum(lam0)
      #extended to all event times
      Lam0 <- rep(Lam0,rutd)
      if(data$stat[1]==0) Lam0 <- c(0,Lam0)
      #extended to all exit times
      Lam0 <- rep(Lam0,rtd)
      # save original object:
      Lam00 <- Lam0
      # for time dependent covariates and left-truncated individuals: replace by the difference
      if(any(start!=0)){
        timehaz <- data.frame(time=sort(data$Y), Lam0_2=Lam0)
        timehaz_tmp <- data.frame(time=unique(data$start), Lam0_2=NA)
        timehaz <- rbind(timehaz, timehaz_tmp)
        timehaz <- timehaz[order(timehaz$time),]
        timehaz$Lam0_2 <- mstateNAfix(timehaz$Lam0_2, 0)
        timehaz <- timehaz[!duplicated(timehaz$time),]

        # Prepare object so that you can calculate Lam0_event_time - Lam0_entry_time
        data_lt <- cbind(data, Lam0, id_0=1:nrow(data))
        data_lt <- merge(data_lt, timehaz,
                         by.x='start', by.y='time', all.x = TRUE)
        data_lt <- data_lt[order(data_lt$id_0),]

        # Edit Lam0:
        Lam0[start!=0] <- data_lt$Lam0[start!=0] - data_lt$Lam0_2[start!=0]
        # Lam0[start!=0] <- Lam0[start!=0] - Lam0[wstart[start!=0]]
      }

      #the smooting of lam0
      if(bwin[1]!=0){
        lam0s <- krn%*%lam0
      } else{
        lam0s <- lam0/difftu
      }

      #extended to all event times
      lam0s <- rep(lam0s,rutd)

      #compute Nie, only for those with unknown hazard
      Nie[cause==2] <- as.vector(lam0s*eb/(a+lam0s*eb))[cause==2]

      #likelihood calculation - manjka ti se likelihood za nicelni model!!!

      lam0 <- rep(lam0,rutd)

      likely <- sum(log(a0 + lam0*eb)) - sum(data$ds + Lam0*ebx)

      if(p>0){
        tempind <- Nie<=0|Nie>=1
        if(any(tempind)){
          if(any(Nie<=0))Nie[Nie<=0] <- tol
          if(any(Nie>=1))Nie[Nie>=1] <- 1-tol
          #if(which(tempind)!=nev)warning("Weights smaller than 0")
          #if(any(is.na( match(which(tempind),c(1,nev)) )))browser()
        }
        if(nded==0) break()
        databig$wei[1:nded] <- c(Nie[cause==2],1-Nie[cause==2])
        bd <- abs(b-b00)
        if(max(bd)< tol) break()
      }
      #early stopping time for no covariates???
    }
  }
  iter <- i
  #if (maxiter > 1& iter>=maxiter)
  #        warning("Ran out of iterations and did not converge")
  if(p>0){
    if(nded!=0){
      resi <- resid(fit,type="schoenfeld")
      if(!is.null(dim(resi)))resi <- resi[1:(nded/2),]
      else resi <- resi[1:(nded/2)]
      swei <- fit$weights[1:(nded/2)]

      if(is.null(dim(resi))) fishem <- sum((resi^2*swei*(1-swei)))
      else {
        fishem <- apply(resi,1,function(x)outer(x,x))
        fishem <- t(t(fishem)*swei*(1-swei))
        fishem <- matrix(apply(fishem,1,sum),ncol=m)
      }
    }

    else fishem <- 0
    fishcox <- solve(fit$var)
    fisher <- fishcox - fishem
    fit$var <- solve(fisher)
    names(fit$coefficients)<-varnames
    fit$lambda0 <- lam0s
  }
  else fit <- list(lambda0 = lam0s)

  # lam0.ns -> lam0u: editan (smoothing) [manj editan k lam0s]
  # lambda0 -> lam0s: editan (smoothing)
  # Lambda0 -> Lam0: cumsum(lam0) - lam0 je enako lam0u

  fit$lambda0 <- fit$lambda0[utimes]
  fit$Lambda0 <- Lam00[udtimes]
  fit$times <- utd
  fit$Nie <- Nie
  fit$bwin <- bwin
  fit$iter <- i
  class(fit) <- c("rsadd",class(fit))
  fit$loglik <- c(likely0,likely)
  fit$lam0.ns <- lam0u
  fit
}



em <- function (rform, init, control, bwin)
{
  data <- rform$data
  n <- nrow(data)
  p <- rform$m
  ord_id <- order(data$Y)
  rform$cause <- rform$cause[ord_id]
  data <- data[ord_id, ]
  fk <- (attributes(rform$ratetable)$factor != 1)
  nfk <- length(fk)
  nev <- length(data$Y[data$stat == 1])
  data$a <- rep(NA, n)
  xx <- exp_prep(data[, 4:(nfk + 3),drop=FALSE], data$Y - data$start, rform$ratetable)
  # The cumulative population hazard of dying at time Y:
  data$ds <- -log(xx)
  data1 <- data
  data1[, 4:(nfk + 3)] <- data[, 4:(nfk + 3)] + data$Y %*% t(fk)
  xx <- exp_prep(data1[data1$stat == 1, 4:(nfk + 3),drop=FALSE], 1, rform$ratetable)
  # The population hazard of dying in the following day (for individuals that had an event):
  data$a[data$stat == 1] <- -log(xx)

  if (p > 0) {
    if (!missing(init) && !is.null(init)) {
      if (length(init) != p)
        stop("Wrong length for inital values")
    }
    else init <- rep(0, p)
    beta <- matrix(init, p, 1)
  }
  pr.time<-proc.time()[3]

  Nie <- rep(.5,sum(data$stat==1))
  Nie[rform$cause[data$stat==1]<2] <-  rform$cause[data$stat==1][rform$cause[data$stat==1]<2]

  #NEW IN 2.05
  varstart <- 3+nfk+1		#first column of covariates
  varstop <- 3+nfk+p		#last column of covariates


  if(missing(bwin))bwin <- -1
  if(bwin<0){

    if(p>0)data1 <- data[,-c(varstart:varstop)]    		#NEW IN 2.05
    else  data1 <- data
    nfk <- length(attributes(rform$ratetable)$dimid)
    names(data)[4:(3+nfk)] <- attributes(rform$ratetable)$dimid
    expe <- rs.surv(Surv(Y,stat)~1,data,ratetable=rform$ratetable,method="ederer2")
    esurv <- -log(expe$surv[expe$n.event!=0])
    if(esurv[length(esurv)]==Inf)esurv[length(esurv)] <-  esurv[length(esurv)-1]
    x <- seq(.1,3,length=5)
    dif <- rep(NA,5)
    options(warn=-1)
    diter <- max(round(max(data$Y)/356.24),3)
    for(it in 1:5){
      fit <- rsfitterem(data1,NULL,diter,rform$ratetable,control$epsilon,x[it],0,rform$cause,Nie)
      dif[it] <- sum((esurv-fit$Lambda0)^2)
    }
    wh <- which.min(dif)
    if(wh==1)x <- seq(x[wh],x[wh+1]-.1,length=5)
    else if(wh==5)x <- c(x, max(data$Y)/ max(diff(data$Y)))
    if(wh!=1)
      x <- seq(x[wh-1]+.1,x[wh+1]-.1,length=5)
    dif <- rep(NA,5)

    for(it in 1:5){
      fit <- rsfitterem(data1,NULL,diter,rform$ratetable,control$epsilon,x[it],0,rform$cause,Nie)
      dif[it] <- sum((esurv-fit$Lambda0)^2)
    }
    options(warn=0)
    Nie <- fit$Nie
    bwin <- x[which.min(dif)]
  }

  # if(any(data$start != 0)) browser()

  fit <- rsfitterem(data, beta, control$maxit, rform$ratetable,
                    control$epsilon, bwin, p, rform$cause,Nie)

  Nie <- rep(0,nrow(data))
  Nie[data$stat==1] <- fit$Nie
  fit$Nie <- Nie[order(ord_id)]
  fit$bwin <- list(bwin=fit$bwin,bwinfac=bwin)
  fit
}




#' Fit an Additive model for Relative Survival
#'
#' The function fits an additive model to the data. The methods implemented are
#' the maximum likelihood method, the semiparametric method, a glm model with a
#' \code{binomial} error and a glm model with a \code{poisson} error.
#'
#' NOTE: The follow-up time must be specified in days. The \code{ratetable}
#' being used may have different variable names and formats than the user's
#' data set, this is dealt with by the \code{rmap} argument. For example, if
#' age is in years in the data set but in days in the \code{ratetable} object,
#' age=age*365.241 should be used. The calendar year can be in any date format
#' (Date and POSIXt are allowed), the date formats in the
#' \code{ratetable} and in the data may differ.
#'
#' The maximum likelihood method and both glm methods assume a fully parametric
#' model with a piecewise constant baseline excess hazard function. The
#' intervals on which the baseline is assumed constant should be passed via
#' argument \code{int}. The EM method is semiparametric, i.e. no assumptions
#' are made for the baseline hazard and therefore no intervals need to be
#' specified.
#'
#' The methods using glm are methods for grouped data. The groups are formed
#' according to the covariate values. This should be taken into account when
#' fitting a model. The glm method returns life tables for groups specified by
#' the covariates in \code{groups}.
#'
#' The EM method output includes the smoothed baseline excess hazard
#' \code{lambda0}, the cumulative baseline excess hazard \code{Lambda0} and
#' \code{times} at which they are estimated. The individual probabilites of
#' dying due to the excess risk are returned as \code{Nie}.  The EM method
#' fitting procedure requires some local smoothing of the baseline excess
#' hazard. The default \code{bwin=-1} value lets the function find an
#' appropriate value for the smoothing band width. While this ensures an
#' unbiased estimate, the procedure time is much longer. As the value found by
#' the function is independent of the covariates in the model, the value can be
#' read from the output (\code{bwinfac}) and used for refitting different
#' models to the same data to save time.
#'
#' @param formula a formula object, with the response as a \code{Surv} object
#' on the left of a \code{~} operator, and, if desired, terms separated by the
#' \code{+} operator on the right. \code{Surv(start,stop,event)} outcomes
#' are also possible for time-dependent covariates and left-truncation for
#' \code{method='EM'}.
#'
#' NOTE: The follow-up time must be in days.
#' @param data a data.frame in which to interpret the variables named in the
#' \code{formula}.
#' @param ratetable a table of event rates, organized as a \code{ratetable}
#' object, such as \code{slopop}.
#' @param int either a single value denoting the number of follow-up years or a
#' vector specifying the intervals (in years) in which the hazard is constant
#' (the times that are bigger than \code{max(int)} are censored. If missing,
#' only one interval (from time 0 to maximum observation time) is assumed.  The
#' EM method does not need the intervals, only the maximum time can be
#' specified (all times are censored after this time point).
#' @param na.action a missing-data filter function, applied to the model.frame,
#' after any subset argument has been used.  Default is
#' \code{options()$na.action}.
#' @param method \code{glm.bin} or \code{glm.poi} for a glm model, \code{EM}
#' for the EM algorithm and \code{max.lik} for the maximum likelihood model
#' (default).
#' @param init vector of initial values of the iteration.  Default initial
#' value is zero for all variables.
#' @param bwin controls the bandwidth used for smoothing in the EM algorithm.
#' The follow-up time is divided into quartiles and \code{bwin} specifies a
#' factor by which the maximum between events time length on each interval is
#' multiplied. The default \code{bwin=-1} lets the function find an appropriate
#' value. If \code{bwin=0}, no smoothing is applied.
#' @param centered if \code{TRUE}, all the variables are centered before
#' fitting and the baseline excess hazard is calculated accordingly. Default is
#' \code{FALSE}.
#' @param cause A vector of the same length as the number of cases. \code{0}
#' for population deaths, \code{1} for disease specific deaths, \code{2}
#' (default) for unknown. Can only be used with the \code{EM} method.
#' @param control a list of parameters for controlling the fitting process.
#' See the documentation for \code{glm.control} for details.
#' @param rmap an optional list to be used if the variables are not organized
#' and named in the same way as in the \code{ratetable} object. See details
#' below.
#' @param ... other arguments will be passed to \code{glm.control}.
#' @return An object of class \code{rsadd}. In the case of
#' \code{method="glm.bin"} and \code{method="glm.poi"} the class also inherits
#' from \code{glm} which inherits from the class \code{lm}. Objects of this
#' class have methods for the functions \code{print} and \code{summary}. An
#' object of class \code{rsadd} is a list containing at least the following
#' components: \item{data}{the data as used in the model, along with the
#' variables defined in the rate table} \item{ratetable}{the ratetable used.}
#' \item{int}{the maximum time (in years) used. All the events at and after
#' this value are censored.} \item{method}{the fitting method that was used.}
#' \item{linear.predictors}{the vector of linear predictors, one per subject.}
#' @seealso \code{\link{rstrans}}, \code{\link{rsmul}}
#' @references Package. Pohar M., Stare J. (2006) "Relative survival analysis
#' in R." Computer Methods and Programs in Biomedicine, \bold{81}: 272--278
#'
#' Relative survival: Pohar, M., Stare, J. (2007) "Making relative survival
#' analysis relatively easy." Computers in biology and medicine, \bold{37}:
#' 1741--1749.
#'
#' EM algorithm: Pohar Perme M., Henderson R., Stare, J. (2009) "An approach to
#' estimation in relative survival regression." Biostatistics, \bold{10}:
#' 136--146.
#' @keywords survival
#' @examples
#'
#' data(slopop)
#' data(rdata)
#' #fit an additive model
#' #note that the variable year is given in days since 01.01.1960 and that
#' #age must be multiplied by 365.241 in order to be expressed in days.
#' fit <- rsadd(Surv(time,cens)~sex+as.factor(agegr)+ratetable(age=age*365.241),
#' 	    ratetable=slopop,data=rdata,int=5)
#'
#' #check the goodness of fit
#' rs.br(fit)
#'
#' #use the EM method and plot the smoothed baseline excess hazard
#' fit <- rsadd(Surv(time,cens)~sex+age,rmap=list(age=age*365.241),
#' 	 ratetable=slopop,data=rdata,int=5,method="EM")
#' sm <- epa(fit)
#' plot(sm$times,sm$lambda,type="l")
#'
rsadd <- function (formula = formula(data), data = parent.frame(), ratetable = relsurv::slopop,
                   int, na.action, method = "max.lik", init, bwin, centered = FALSE,
                   cause, control, rmap, ...) # rmap.subs=TRUE,
{
  call <- match.call()
  if (missing(control))
    control <- glm.control(...)

  if (!missing(rmap)) {
    rmap_tmp <- substitute(rmap)
    if(inherits(rmap_tmp, "call")){
      rmap <- rmap_tmp
    }
  }

  if(!missing(cause)){								#NEW: ce cause ne manjka, ga preverim in dodam kot spremenljivko
    if (length(cause) != nrow(data))
      stop("Length of cause does not match data dimensions")
    data$cause <- cause
    rform <- rformulate(formula, data, ratetable, na.action, rmap, int, centered, cause)
  }
  else{ #no cause
    rform <- rformulate(formula,data, ratetable, na.action, rmap, int, centered)
  }

  if (method == "EM") {
    if (!missing(int)) {
      if (length(int) > 1 | any(int <= 0))
        stop("Invalid value of 'int'")
    }
  }
  else {
    if (missing(int))
      int <- c(0,ceiling(max(rform$Y/365.241)))
    if (length(int) == 1) {
      if (int <= 0)
        stop("The value of 'int' must be positive ")
      int <- 0:int
    }
    else if (int[1] != 0)
      stop("The first interval in 'int' must start with 0")
  }
  method <-  match.arg(method,c("glm.bin","glm.poi","max.lik","EM"))

  if (method == "glm.bin" | method == "glm.poi")
    fit <- glmxp(rform = rform, interval = int, method = method,
                 control = control)
  else if (method == "max.lik")
    fit <- maxlik(rform = rform, interval = int, init = init,
                  control = control)
  else if (method == "EM")
    fit <- em(rform, init, control, bwin)
  fit$call <- call
  fit$formula <- formula
  fit$data <- rform$data
  fit$ratetable <- rform$ratetable
  fit$n <- nrow(rform$data)
  if (length(rform$na.action))
    fit$na.action <- rform$na.action
  fit$y <- rform$Y.surv
  fit$method <- method
  if (method == "EM") {
    if (!missing(int))
      fit$int <- int
    else fit$int <- ceiling(max(rform$Y[rform$status == 1])/365.241)
    fit$terms <- rform$Terms
    if(centered)fit$mvalue <- rform$mvalue
  }
  if (method == "max.lik") {
    fit$terms <- rform$Terms
  }
  if (rform$m > 0)
    fit$linear.predictors <- as.matrix(rform$X) %*% fit$coef[1:ncol(rform$X)]

  if (!missing(rmap)) {
    fit$rmap <- rmap
  }
  fit
}


maxlik <- function (rform, interval, subset, init, control)
{
  data <- rform$data
  max.time <- max(data$Y)/365.241
  if (max.time < max(interval))
    interval <- interval[1:(sum(max.time > interval) + 1)]
  fk <- (attributes(rform$ratetable)$factor != 1)
  nfk <- length(fk)
  data <- cbind(data, offset = rform$offset)
  data <- survsplit(data, cut = interval[-1] * 365.241, end = "Y",
                    event = "stat", start = "start", episode = "epi", interval = interval)
  del <- which(data$start==data$Y)
  if(length(del))    data <- data[-del,]
  offset <- data$offset
  data$offset <- NULL
  d.int <- diff(interval)
  data[, 4:(nfk + 3)] <- data[, 4:(nfk + 3)] + data$start %*%
    t(fk)
  data$lambda <- rep(0, nrow(data))
  nsk <- nrow(data[data$stat == 1, ])
  xx <- exp_prep(data[data$stat == 1, 4:(nfk + 3),drop=FALSE] + (data[data$stat ==
                                                                        1, ]$Y - data[data$stat == 1, ]$start) %*% t(fk), 1,  rform$ratetable)
  data$lambda[data$stat == 1] <- -log(xx) * 365.241
  xx <- exp_prep(data[, 4:(nfk + 3),drop=FALSE], data$Y - data$start, rform$ratetable)
  data$epi <- NULL
  data$ds <- -log(xx)
  data$Y <- data$Y/365.241
  data$start <- data$start/365.241
  data <- data[, -(4:(3 + nfk))]
  intn <- length(interval[-1])
  m <- rform$m
  p <-  m + intn
  if (!missing(init) && !is.null(init)) {
    if (length(init) != p)
      stop("Wrong length for inital values")
  }
  else init <- rep(0, p)
  if(m>0){
    init0 <- init[-(1:m)]
    data1 <- data[,-(4:(3+m))]
  }
  else{
    init0 <- init
    data1 <- data
  }
  fit0 <- lik.fit(data1, 0, intn, init0,  control, offset)
  if(m>0){
    init[-(1:m)] <- fit0$coef
    fit <- lik.fit(data, m, intn, init,  control, offset)
  }
  else fit <- fit0
  fit$int <- interval
  class(fit) <- "rsadd"
  fit$times <- fit$int*365.241						#dodano za potrebe rs.surv.rsadd
  fit$Lambda0 <- cumsum(c(0, exp(fit$coef[(m+1):p])*diff(fit$int)  ))
  fit
}

lik.fit <- function (data, m, intn, init, control, offset)
{
  n <- dim(data)[1]
  varpos <- 4:(3 + m + intn)
  x <- data[, varpos]
  varnames <- names(data)[varpos]
  lbs <- names(x)
  x <- as.matrix(x)
  p <- length(varpos)
  d <- data$stat
  ds <- data$ds
  h <- data$lambda
  y <- data$Y - data$start
  maxiter <- control$maxit
  if (!missing(init) && !is.null(init)) {
    if (length(init) != p)
      stop("Wrong length for inital values")
  }
  else init <- rep(0, p)
  b <- matrix(init, p, 1)
  b0 <- b
  fit <- mlfit(b, p, x, offset, d, h, ds, y, maxiter, control$epsilon)
  if (maxiter > 1 & fit$nit >= maxiter) {
    values <- apply(data[data$stat==1,varpos,drop=FALSE],2,sum)			#NEW: deluje tudi, ce je ratetable eno-dimenzionalen
    problem <- which.min(values)
    outmes <- "Ran out of iterations and did not converge"
    if(values[problem]==0)tzero <- ""
    else tzero <- "only "
    if(values[problem]<5){
      if(!is.na(strsplit(names(values)[problem],"fu")[[1]][2]))outmes <- paste(outmes, "\n This may be due to the fact that there are ",tzero, values[problem], " events on interval",strsplit(names(values)[problem],"fu")[[1]][2],"\n You can use the 'int' argument to change the follow-up intervals in which the baseline excess hazard is assumed constant",sep="")
      else outmes <- paste(outmes, "\n This may be due to the fact that there are ",tzero, values[problem], " events for covariate value ",names(values)[problem],sep="")
    }
    warning(outmes)
  }
  b <- as.vector(fit$b)
  names(b) <- varnames
  fit <- list(coefficients = b, var = -solve(fit$sd), iter = fit$nit,
              loglik = fit$loglik)
  fit
}





#' Split a Survival Data Set at Specified Times
#'
#' Given a survival data set and a set of specified cut times, the function
#' splits each record into multiple records at each cut time. The new data set
#' is be in \code{counting process} format, with a start time, stop time, and
#' event status for each record. More general than \code{survSplit} as it also
#' works with the data already in the \code{counting process} format.
#'
#'
#' @param data data frame.
#' @param cut vector of timepoints to cut at.
#' @param end character string with name of event time variable.
#' @param event character string with name of censoring indicator.
#' @param start character string with name of start variable (will be created
#' if it does not exist).
#' @param id character string with name of new id variable to create
#' (optional).
#' @param zero If \code{start} doesn't already exist, this is the time that the
#' original records start. May be a vector or single value.
#' @param episode character string with name of new episode variable
#' (optional).
#' @param interval this argument is used by \code{max.lik} function
#' @return New, longer, data frame.
#' @seealso \code{\link[survival:survSplit]{survival::survSplit}}.
#' @keywords survival
survsplit <- function (data, cut, end, event, start, id = NULL, zero = 0,
                       episode = NULL, interval = NULL)
{
  ntimes <- length(cut)
  n <- nrow(data)
  p <- ncol(data)
  if (length(interval) > 0) {
    ntimes <- ntimes - 1
    sttime <- c(rep(0, n), rep(cut[-length(cut)], each = n))
    endtime <- rep(cut, each = n)
  }
  else {
    endtime <- rep(c(cut, Inf), each = n)
    sttime <- c(rep(0, n), rep(cut, each = n))
  }
  newdata <- lapply(data, rep, ntimes + 1)
  eventtime <- newdata[[end]]
  if (start %in% names(data))
    starttime <- newdata[[start]]
  else starttime <- rep(zero, length = (ntimes + 1) * n)
  starttime <- pmax(sttime, starttime)
  epi <- rep(0:ntimes, each = n)
  if (length(interval) > 0)
    status <- ifelse(eventtime <= endtime & eventtime >=
                       starttime, newdata[[event]], 0)
  else status <- ifelse(eventtime <= endtime & eventtime >
                          starttime, newdata[[event]], 0)
  endtime <- pmin(endtime, eventtime)
  if (length(interval) > 0)
    drop <- (starttime > endtime) | (starttime == endtime &
                                       status == 0)
  else drop <- starttime >= endtime
  newdata <- do.call("data.frame", newdata)
  newdata <- newdata[!drop, ]
  newdata[, start] <- starttime[!drop]
  newdata[, end] <- endtime[!drop]
  newdata[, event] <- status[!drop]
  if (!is.null(id))
    newdata[, id] <- rep(rownames(data), ntimes + 1)[!drop]
  fu <- NULL
  if (length(interval) > 2) {
    for (it in 1:length(interval[-1])) {
      drop1 <- sum(!drop[1:(it * n - n)])
      drop2 <- sum(!drop[(it * n - n + 1):(it * n)])
      drop3 <- sum(!drop[(it * n + 1):(length(interval[-1]) *
                                         n)])
      if (it == 1)
        fu <- cbind(fu, c(rep(1, drop2), rep(0, drop3)))
      else if (it == length(interval[-1]))
        fu <- cbind(fu, c(rep(0, drop1), rep(1, drop2)))
      else fu <- cbind(fu, c(rep(0, drop1), rep(1, drop2),
                             rep(0, drop3)))
    }
    fu <- as.data.frame(fu)
    names(fu) <- c(paste("fu [", interval[-length(interval)],
                         ",", interval[-1], ")", sep = ""))
    newdata <- cbind(newdata, fu)
  }
  else if (length(interval) == 2) {
    fu <- rep(1, sum(!drop))
    newdata <- cbind(newdata, fu)
    names(newdata)[ncol(newdata)] <- paste("fu [", interval[1],
                                           ",", interval[2], "]", sep = "")
  }
  if (!is.null(episode))
    newdata[, episode] <- epi[!drop]
  newdata
}


glmxp <- function (rform, data, interval, method, control)
{
  if (rform$m == 1)
    g <- as.integer(as.factor(rform$X[[1]]))
  else if (rform$m > 1) {
    gvar <- NULL
    for (i in 1:rform$m) {
      gvar <- append(gvar, rform$X[i])
    }
    tabgr <- as.data.frame(table(gvar))
    tabgr <- tabgr[, 1:rform$m]
    n.groups <- dim(tabgr)[1]
    mat <- do.call("data.frame", gvar)
    names(mat) <- names(tabgr)
    tabgr <- cbind(tabgr, g = as.numeric(row.names(tabgr)))
    mat <- cbind(mat, id = 1:rform$n)
    c <- merge(tabgr, mat)
    g <- c[order(c$id), rform$m + 1]
  }
  else g <- rep(1, rform$n)
  vg <- function(X) {
    n <- dim(X)[1]
    w <- sum((X$event == 0) & (X$fin == 1) & (X$y != 1))
    nd <- sum((X$event == 1) & (X$fin == 1))
    ps <- exp_prep(X[, 4:(nfk + 3),drop=FALSE], t.int, rform$ratetable)
    ld <- n - w/2
    lny <- log(sum(X$y))
    k <- t.int/365.241
    dstar <- sum(-log(ps)/k * X$y)
    ps <- mean(ps)
    if (rform$m == 0)
      data.rest <- X[1, 7 + nfk + rform$m, drop = FALSE]
    else data.rest <- X[1, c((3 + nfk + 1):(3 + nfk + rform$m),
                             7 + nfk + rform$m)]
    cbind(nd = nd, ld = ld, ps = ps, lny = lny, dstar = dstar,
          k = k, data.rest)
  }
  nint <- length(interval)
  if (nint < 2)
    stop("Illegal interval value")
  meje <- interval
  my.fun <- function(x) {
    if (x > 1) {
      x.t <- rep(1, floor(x))
      if (x - floor(x) > 0)
        x.t <- c(x.t, x - floor(x))
      x.t
    }
    else x
  }
  int <- apply(matrix(diff(interval), ncol = 1), 1, my.fun)
  if (is.list(int))
    int <- c(0, cumsum(do.call("c", int)))
  else int <- c(0, cumsum(int))
  int <- int * 365.241
  nint <- length(int)
  X <- cbind(rform$data, grupa = g)
  fk <- (attributes(rform$ratetable)$factor != 1)
  nfk <- length(fk)
  Z <- X[X$start >= int[2], ]
  nz <- dim(Z)[1]
  Z$fin <- rep(0, nz)
  Z$event <- rep(0, nz)
  Z$fu <- rep(0, nz)
  Z$y <- rep(0, nz)
  Z$origstart <- Z$start
  Z$xind <- rep(0, nz)
  if (nrow(Z) > 0)
    Z[, 4:(nfk + 3)] <- Z[, 4:(nfk + 3)] + matrix(Z$start,
                                                  ncol = nfk, byrow = FALSE, nrow = nrow(Z)) * matrix(fk,
                                                                                                      ncol = nfk, byrow = TRUE, nrow = nrow(Z))
  X <- X[X$start < int[2], ]
  X$fin <- (X$Y <= int[2])
  X$event <- X$fin * X$stat
  ford <- eval(substitute(paste("[", a, ",", b, "]", sep = ""),
                          list(a = meje[1], b = meje[2])))
  X$fu <- rep(ford, rform$n - nz)
  t.int <- int[2] - int[1]
  X$y <- (pmin(X$Y, int[2]) - X$start)/365.241
  X$origstart <- X$start
  X$xind <- rep(1, nrow(X))
  gr1 <- by(X, X$grupa, vg)
  grm1 <- do.call("rbind", gr1)
  X <- X[X$fin == 0, ]
  X$start <- rep(int[2], dim(X)[1])
  X <- rbind(X, Z[Z$start < int[3], ])
  Z <- Z[Z$start >= int[3], ]
  temp <- 0
  if (nint > 2) {
    for (i in 3:nint) {
      ni <- dim(X)[1]
      if (ni == 0) {
        temp <- 1
        break
      }
      X$fin <- X$Y <= int[i]
      X$event <- X$fin * X$stat
      l <- sum(int[i - 1] >= meje * 365.241)
      if(l==1)
        ftemp <- eval(substitute(paste("[", a, ",", b, "]", sep = ""),
                                 list(a = meje[l], b = meje[l + 1])))

      else
        ftemp <- eval(substitute(paste("(", a, ",", b, "]", sep = ""),
                                 list(a = meje[l], b = meje[l + 1])))
      ford <- c(ford, ftemp)
      X$fu <- rep(ford[i - 1], ni)
      t.int <- int[i] - int[i - 1]
      index <- X$origstart < int[i - 1]
      index1 <- as.logical(X$xind)
      if (sum(index) > 0)
        X[index, 4:(nfk + 3)] <- X[index, 4:(nfk + 3)] +
        matrix(fk * t.int, ncol = nfk, byrow = TRUE,
               nrow = sum(index))
      X$xind <- rep(1, nrow(X))
      X$y <- (pmin(X$Y, int[i]) - X$start)/365.241
      gr1 <- by(X, X$grupa, vg)
      grm1 <- rbind(grm1, do.call("rbind", gr1))
      X <- X[X$fin == 0, ]
      X$start <- rep(int[i], dim(X)[1])
      if (i == nint)
        break
      X <- rbind(X, Z[Z$start < int[i + 1], ])
      X <- X[X$start != X$Y, ]
      Z <- Z[Z$start >= int[i + 1], ]
    }
    l <- sum(int[i - temp] > meje * 365.241)
    interval <- meje[1:(l + 1)]
  }
  else interval <- meje[1:2]
  grm1$fu <- factor(grm1$fu, levels = unique(ford))
  if (method == "glm.bin") {
    ht <- binomial(link = cloglog)
    ht$link <- "Hakulinen-Tenkanen relative survival model"
    ht$linkfun <- function(mu) log(-log((1 - mu)/ps))
    ht$linkinv <- function(eta) 1 - exp(-exp(eta)) * ps
    ht$mu.eta <- function(eta) exp(eta) * exp(-exp(eta)) *
      ps
    .ps <- ps <- grm1$ps
    #assign(".ps", grm1$ps, envir = .GlobalEnv)
    # ht$initialize <- expression({
    #     n <- y[, 1] + y[, 2]
    #     y <- ifelse(n == 0, 0, y[, 1]/n)
    #     weights <- weights * n
    #     mustart <- (n * y + 0.01)/(n + 0.02)
    #     mustart[(1 - mustart)/data$ps >= 1] <- data$ps[(1 - mustart)/data$ps >=
    #         1] * 0.9
    # })
    if (any(grm1$ld - grm1$nd > grm1$ps * grm1$ld)) {
      n <- sum(grm1$ld - grm1$nd > grm1$ps * grm1$ld)
      g <- dim(grm1)[1]
      warnme <- paste("Observed number of deaths is smaller than the expected in ",
                      n, "/", g, " groups of patients", sep = "")
    }
    else warnme <- ""

    if (length(interval) == 2 & rform$m == 0)
      stop("No groups can be formed")
    if (length(interval) == 1 | length(table(grm1$fu)) ==
        1)
      grm1$fu <- as.integer(grm1$fu)

    y <- ifelse(grm1$ld == 0, 0, grm1$nd/grm1$ld)
    #weights <- weights * grm1$ld
    mustart <- (grm1$ld * y + 0.01)/(grm1$ld + 0.02)
    mustart[(1 - mustart)/grm1$ps >= 1] <- grm1$ps[(1 - mustart)/grm1$ps >=
                                                     1] * 0.9

    if (!length(rform$X))
      local.ht <- glm(cbind(nd, ld - nd) ~ -1 + fu + offset(log(k)),
                      data = grm1, family = ht,mustart=mustart)
    else {
      xmat <- as.matrix(grm1[, 7:(ncol(grm1) - 1)])

      local.ht <- glm(cbind(nd, ld - nd) ~ -1 + xmat +
                        fu + offset(log(k)), data = grm1, family = ht,mustart=mustart)
    }
    names(local.ht[[1]]) <- c(names(rform$X), paste("fu",
                                                    levels(grm1$fu)))
  }
  else if (method == "glm.poi") {
    pot <- poisson()
    pot$link <- "glm relative survival model with Poisson error"
    pot$linkfun <- function(mu) log(mu - dstar)
    pot$linkinv <- function(eta) dstar + exp(eta)
    #assign(".dstar", grm1$dstar, envir = .GlobalEnv)
    if (any(grm1$nd - grm1$dstar < 0)) {
      pot$initialize <- expression({
        if (any(y < 0)) stop(paste("Negative values not allowed for",
                                   "the Poisson family"))
        n <- rep.int(1, nobs)
        #mustart <- pmax(y, .dstar) + 0.1
      })
    }
    if (any(grm1$nd - grm1$dstar < 0)) {
      n <- sum(grm1$nd - grm1$dstar < 0)
      g <- dim(grm1)[1]
      warnme <- paste("Observed number of deaths is smaller than the expected in ",
                      n, "/", g, " groups of patients", sep = "")
    }
    else warnme <- ""
    dstar <- grm1$dstar
    if (length(interval) == 2 & rform$m == 0)
      stop("No groups can be formed")
    if (length(interval) == 1 | length(table(grm1$fu)) ==
        1)
      grm1$fu <- as.integer(grm1$fu)

    mustart <- pmax(grm1$nd, grm1$dstar) + 0.1
    if (!length(rform$X))
      local.ht <- glm(nd ~ -1 + fu, data = grm1, family = pot,
                      offset = grm1$lny,mustart=mustart)
    else {
      xmat <- as.matrix(grm1[, 7:(ncol(grm1) - 1)])
      local.ht <- glm(nd ~ -1 + xmat + fu, data = grm1,
                      family = pot, offset = grm1$lny,mustart=mustart)
    }
    names(local.ht[[1]]) <- c(names(rform$X), paste("fu",
                                                    levels(grm1$fu)))
  }
  else stop(paste("Method '", method, "' not a valid method",
                  sep = ""))
  class(local.ht) <- c("rsadd", class(local.ht))
  local.ht$warnme <- warnme
  local.ht$int <- interval
  local.ht$groups <- local.ht$data
  return(local.ht)
}



#' Calculate Residuals for a "rsadd" Fit
#'
#' Calculates partial residuals for an additive relative survival model.
#'
#'
#' @param object an object inheriting from class \code{rsadd}, representing a
#' fitted additive relative survival model. Typically this is the output from
#' the \code{rsadd} function.
#' @param type character string indicating the type of residual desired.
#' Currently only Schoenfeld residuals are implemented.
#' @param ... other arguments.
#' @return A list of the following values is returned: \item{res}{a matrix
#' containing the residuals for each variable.} \item{varr}{the variance for
#' each residual} \item{varr1}{the sum of \code{varr}.} \item{kvarr}{the
#' derivative of each residual, to be used in \code{rs.zph} function.}
#' \item{kvarr1}{the sum of \code{kvarr}.}
#' @seealso \code{\link{rsadd}}.
#' @references Package. Pohar M., Stare J. (2006) "Relative survival analysis
#' in R." Computer Methods and Programs in Biomedicine, \bold{81}: 272--278
#'
#' Relative survival: Pohar, M., Stare, J. (2007) "Making relative survival
#' analysis relatively easy." Computers in biology and medicine, \bold{37}:
#' 1741--1749.
#'
#' Goodness of fit: Stare J.,Pohar Perme M., Henderson R. (2005) "Goodness of
#' fit of relative survival models." Statistics in Medicine, \bold{24}:
#' 3911--3925.
#' @keywords survival
#' @examples
#'
#' data(slopop)
#' data(rdata)
#' fit <- rsadd(Surv(time,cens)~sex,rmap=list(age=age*365.241),
#'        ratetable=slopop,data=rdata,int=5)
#' sresid <- residuals.rsadd(fit)
#'
residuals.rsadd <- function (object, type = "schoenfeld", ...)
{
  data <- object$data[order(object$data$Y), ]
  ratetable <- object$ratetable
  beta <- object$coef
  start <- data[, 1]
  stop <- data[, 2]
  event <- data[, 3]
  fk <- (attributes(ratetable)$factor != 1)
  nfk <- length(fk)
  n <- nrow(data)
  scale <- 1
  if (object$method == "EM")
    scale <- 365.241
  m <- ncol(data)
  rem <- m - nfk - 3
  interval <- object$int
  int <- ceiling(max(interval))
  R <- data[, 4:(nfk + 3)]
  lp <- matrix(-log(exp_prep(as.matrix(R), 365.241, object$ratetable))/scale, ncol = 1)
  fu <- NULL
  if (object$method == "EM") {
    death.time <- stop[event == 1]
    for (it in 1:int) {
      fu <- as.data.frame(cbind(fu, as.numeric(death.time/365.241 <
                                                 it & (death.time/365.241) >= (it - 1))))
    }
    if(length(death.time)!=length(unique(death.time))){
      utimes <- which(c(1,diff(death.time))!=0)
      razteg <- function(x){
        # x is a 0/1 vector, the output is a vector of length sum(x), with the corresponding rep numbers
        n <- length(x)
        repu <- rep(1,n)
        repu[x==1] <- 0
        repu <- rev(cumsum(rev(repu)))
        repu <- repu[x==1]
        repu <- -diff(c(repu,0))+1
        if(sum(repu)!=n)repu <- c(n-sum(repu),repu)  #ce je prvi cas censoring, bo treba se kej narest??
        repu
      }
      rutd <- rep(0,length(death.time))
      rutd[utimes] <- 1
      rutd <- razteg(rutd)				#from  unique event times to event times
    }
    else rutd <- rep(1,length(death.time))
    lambda0 <- rep(object$lambda0,rutd)
  }
  else {
    pon <- NULL
    for (i in 1:(length(interval) - 1)) {
      width <- ceiling(interval[i + 1]) - floor(interval[i])
      lo <- interval[i]
      hi <- min(interval[i + 1], floor(interval[i]) + 1)
      for (j in 1:width) {
        fu <- as.data.frame(cbind(fu, as.numeric(stop/365.241 <
                                                   hi & stop/365.241 >= lo)))
        names(fu)[ncol(fu)] <- paste("fu", lo, "-", hi,
                                     sep = "")
        if (j == width) {
          pon <- c(pon, sum(fu[event == 1, (ncol(fu) -
                                              width + 1):ncol(fu)]))
          break()
        }
        else {
          lo <- hi
          hi <- min(interval[i + 1], floor(interval[i]) +
                      1 + j)
        }
      }
    }
    m <- ncol(data)
    data <- cbind(data, fu)
    rem <- m - nfk - 3
    lambda0 <- rep(exp(beta[rem + 1:(length(interval) - 1)]),
                   pon)
    fu <- fu[event == 1, , drop = FALSE]
    beta <- beta[1:rem]
  }
  if (int >= 2) {
    for (j in 2:int) {
      R <- R + matrix(fk * 365.241, ncol = ncol(R), byrow = TRUE,
                      nrow = n)
      xx <- exp_prep(R, 365.241, object$ratetable)
      lp <- cbind(lp, -log(xx)/scale)
    }
  }
  z <- as.matrix(data[, (4 + nfk):m])
  out <- resid.com(start, stop, event, z, beta, lp, lambda0,
                   fu, n, rem, int, type)
  out
}

resid.com <- function (start, stop, event, z, beta, lp, lambda0, fup, n, rem,
                       int, type)
{
  le <- exp(z %*% beta)
  olp <- if (int > 1)
    apply(lp[n:1, ], 2, cumsum)[n:1, ]
  else matrix(cumsum(lp[n:1])[n:1], ncol = 1)
  ole <- cumsum(le[n:1])[n:1]
  lp.st <- lp[order(start), , drop = FALSE]
  le.st <- le[order(start), , drop = FALSE]
  starter <- sort(start)
  starter1 <- c(starter[1], starter[-length(starter)])
  index <- c(TRUE, (starter != starter1)[-1])
  starter <- starter[index]
  val1 <- apply(matrix(starter, ncol = 1), 1, function(x, Y) sum(x >=
                                                                   Y), stop)
  val1 <- c(val1[1], diff(val1), length(stop) - val1[length(val1)])
  olp.st <- (apply(lp.st[n:1, , drop = FALSE], 2, cumsum)[n:1,
                                                          , drop = FALSE])[index, , drop = FALSE]
  olp.st <- apply(olp.st, 2, function(x) rep(c(x, 0), val1))
  olp <- olp - olp.st
  olp <- olp[event == 1, ]
  olp <- apply(fup * olp, 1, sum)
  ole.st <- cumsum(le.st[n:1])[n:1][index]
  ole.st <- rep(c(ole.st, 0), val1)
  ole <- ole - ole.st
  ole <- ole[event == 1] * lambda0
  s0 <- ole + olp
  sc <- NULL
  zb <- NULL
  kzb <- NULL
  f1 <- function(x) rep(mean(x), length(x))
  f2 <- function(x) apply(x, 2, f1)
  f3 <- function(x) apply(x, 1:2, f1)
  ties <- length(unique(stop[event == 1])) != length(stop[event ==
                                                            1])
  for (k in 1:rem) {
    zlp <- apply((z[, k] * lp)[n:1, , drop = FALSE], 2, cumsum)[n:1,
                                                                , drop = FALSE]
    zlp.st <- (apply((z[, k] * lp.st)[n:1, , drop = FALSE],
                     2, cumsum)[n:1, , drop = FALSE])[index, , drop = FALSE]
    zlp.st <- apply(zlp.st, 2, function(x) rep(c(x, 0), val1))
    zlp <- zlp - zlp.st
    zlp <- zlp[event == 1, , drop = FALSE]
    zlp <- apply(fup * zlp, 1, sum)
    zle <- cumsum((z[, k] * le)[n:1])[n:1]
    zle.st <- cumsum((z[, k] * le.st)[n:1])[n:1][index]
    zle.st <- rep(c(zle.st, 0), val1)
    zle <- zle - zle.st
    zle <- zle[event == 1]
    zle <- zle * lambda0
    s1 <- zle + zlp
    zb <- cbind(zb, s1/s0)
    kzb <- cbind(kzb, zle/s0)
  }
  s1ties <- cbind(zb, kzb)
  if (ties) {
    s1ties <- by(s1ties, stop[event == 1], f2)
    s1ties <- do.call("rbind", s1ties)
  }
  zb <- s1ties[, 1:rem, drop = FALSE]
  kzb <- s1ties[, -(1:rem), drop = FALSE]
  sc <- z[event == 1, , drop = FALSE] - zb
  row.names(sc) <- stop[event == 1]
  out.temp <- function(x) outer(x, x, FUN = "*")
  krez <- rez <- array(matrix(NA, ncol = rem, nrow = rem),
                       dim = c(rem, rem, sum(event == 1)))
  for (a in 1:rem) {
    for (b in a:rem) {
      zzlp <- apply((z[, a] * z[, b] * lp)[n:1, , drop = FALSE],
                    2, cumsum)[n:1, , drop = FALSE]
      zzlp.st <- (apply((z[, a] * z[, b] * lp.st)[n:1,
                                                  , drop = FALSE], 2, cumsum)[n:1, , drop = FALSE])[index,
                                                                                                    , drop = FALSE]
      zzlp.st <- apply(zzlp.st, 2, function(x) rep(c(x,
                                                     0), val1))
      zzlp <- zzlp - zzlp.st
      zzlp <- zzlp[event == 1, , drop = FALSE]
      zzlp <- apply(fup * zzlp, 1, sum)
      zzle <- cumsum((z[, a] * z[, b] * le)[n:1])[n:1]
      zzle.st <- cumsum((z[, a] * z[, b] * le.st)[n:1])[n:1][index]
      zzle.st <- rep(c(zzle.st, 0), val1)
      zzle <- zzle - zzle.st
      zzle <- zzle[event == 1]
      zzle <- zzle * lambda0
      s2 <- zzlp + zzle
      s20 <- s2/s0
      ks20 <- zzle/s0
      s2ties <- cbind(s20, ks20)
      if (ties) {
        s2ties <- by(s2ties, stop[event == 1], f2)
        s2ties <- do.call("rbind", s2ties)
      }
      rez[a, b, ] <- rez[b, a, ] <- s2ties[, 1]
      krez[a, b, ] <- krez[b, a, ] <- s2ties[, 2]
    }
  }
  juhu <- apply(zb, 1, out.temp)
  if (is.null(dim(juhu)))
    juhu1 <- array(data = matrix(juhu, ncol = a), dim = c(a,
                                                          a, length(zb[, 1])))
  else juhu1 <- array(data = apply(juhu, 2, matrix, ncol = a),
                      dim = c(a, a, length(zb[, 1])))
  varr <- rez - juhu1
  kjuhu <- apply(cbind(zb, kzb), 1, function(x) outer(x[1:rem],
                                                      x[-(1:rem)], FUN = "*"))
  if (is.null(dim(kjuhu)))
    kjuhu1 <- array(data = matrix(kjuhu, ncol = rem), dim = c(rem,
                                                              rem, length(zb[, 1])))
  else kjuhu1 <- array(data = apply(kjuhu, 2, matrix, ncol = rem),
                       dim = c(rem, rem, length(zb[, 1])))
  kvarr <- krez - kjuhu1
  for (i in 1:dim(varr)[1]) varr[i, i, which(varr[i, i, ] <
                                               0)] <- 0
  for (i in 1:dim(kvarr)[1]) kvarr[i, i, which(kvarr[i, i,
  ] < 0)] <- 0
  varr1 <- apply(varr, 1:2, sum)
  kvarr1 <- apply(kvarr, 1:2, sum)
  if (type == "schoenfeld")
    out <- list(res = sc, varr1 = varr1, varr = varr, kvarr = kvarr,
                kvarr1 = kvarr1)
  out
}






#' Test the Proportional Hazards Assumption for Relative Survival Regression
#' Models
#'
#' Test the proportional hazards assumption for relative survival models
#' (\code{rsadd}, \code{rsmul} or \code{rstrans}) by forming a Brownian Bridge.
#'
#'
#' @aliases rs.br plot.rs.br print.rs.br
#' @param fit the result of fitting a relative survival model, using the
#' \code{rsadd}, \code{rsmul} or \code{rstrans} function.
#' @param sc partial residuals calculated by the \code{resid} function. This is
#' used to save time if several tests are to be calculated on these residuals
#' and can otherwise be omitted.
#' @param rho a number controlling the weigths of residuals. The weights are
#' the number of individuals at risk at each event time to the power
#' \code{rho}. The default is \code{rho=0}, which sets all weigths to 1.
#' @param test a character string specifying the test to be performed on
#' Brownian bridge. Possible values are \code{"max"} (default), which tests the
#' maximum absolute value of the bridge, and \code{cvm}, which calculates the
#' Cramer Von Mises statistic.
#' @param global should a global Brownian bridge test be performed, in addition
#' to the per-variable tests
#' @return an object of class \code{rs.br}. This function would usually be
#' followed by both a print and a plot of the result. The plot gives a Brownian
#' bridge for each of the variables. The horizontal lines are the 95% and 99%
#' confidence intervals for the maximum absolute value of the Brownian bridge
#' @seealso \code{\link{rsadd}}, \code{rsmul}, \code{rstrans},
#' \code{\link{resid}}.
#' @references Goodness of fit: Stare J.,Pohar Perme M., Henderson R. (2005)
#' "Goodness of fit of relative survival models." Statistics in Medicine,
#' \bold{24}: 3911--3925.
#'
#' Package. Pohar M., Stare J. (2006) "Relative survival analysis in R."
#' Computer Methods and Programs in Biomedicine, \bold{81}: 272--278
#'
#' Relative survival: Pohar, M., Stare, J. (2007) "Making relative survival
#' analysis relatively easy." Computers in biology and medicine, \bold{37}:
#' 1741--1749.
#' @keywords survival
#' @examples
#'
#' data(slopop)
#' data(rdata)
#' fit <- rsadd(Surv(time,cens)~sex,rmap=list(age=age*365.241),
#' 		ratetable=slopop,data=rdata,int=5)
#' rsbr <- rs.br(fit)
#' rsbr
#' plot(rsbr)
#'
rs.br <- function (fit, sc, rho = 0, test = "max", global = TRUE)
{
  test <- match.arg(test,c("max","cvm"))
  if (inherits(fit, "rsadd")) {
    if (missing(sc))
      sc <- resid(fit, "schoenfeld")
    sresid <- sc$res
    varr <- sc$varr
    sresid <- as.matrix(sresid)
  }
  else {
    coef <- fit$coef
    options(warn = -1)
    sc <- coxph.detail(fit)
    options(warn = 0)
    sresid <- sc$score
    varr <- sc$imat
    if (is.null(dim(varr)))
      varr <- array(varr, dim = c(1, 1, length(varr)))
    sresid <- as.matrix(sresid)
  }
  if (inherits(fit, "coxph")) {
    if(is.null(fit$data)){
      temp <- fit$y
      class(temp) <- "matrix"
      if(ncol(fit$y)==2)temp <- data.frame(rep(0,nrow(fit$y)),temp)
      if(is.null(fit$x))stop("The coxph model should be called with x=TRUE argument")
      fit$data <- data.frame(temp,fit$x)
      names(fit$data)[1:3] <- c("start","Y","stat")
    }
  }
  data <- fit$data[order(fit$data$Y), ]
  time <- data$Y[data$stat == 1]
  ties <- (length(unique(time)) != length(time))
  keep <- 1:(ncol(sresid))
  options(warn = -1)
  scaled <- NULL
  varnova <- NULL
  if (ncol(sresid) == 1) {
    varr <- varr[1, 1, ]
    scaled <- sresid/sqrt(varr)
  }
  else { for (i in 1:ncol(sresid)) varnova <- cbind(varnova,varr[i,i,])
  scaled <- sresid/sqrt(varnova)
  }

  options(warn = 0)
  nvar <- ncol(sresid)
  survfit <- getFromNamespace("survfit", "survival")
  temp <- survfit(fit$y~1, type = "kaplan-meier")
  n.risk <- temp$n.risk
  n.time <- temp$time
  if (temp$type == "right") {
    cji <- matrix(fit$y, ncol = 2)
    n.risk <- n.risk[match(cji[cji[, 2] == 1, 1], n.time)]
  }
  else {
    cji <- matrix(fit$y, ncol = 3)
    n.risk <- n.risk[match(cji[cji[, 3] == 1, 2], n.time)]
  }
  n.risk <- sort(n.risk, decreasing = TRUE)
  varnames <- names(fit$coef)[keep]
  u2 <- function(bb) {
    n <- length(bb)
    1/n * (sum(bb^2) - sum(bb)^2/n)
  }
  wc <- function(x, k = 1000) {
    a <- 1
    for (i in 1:k) a <- a + 2 * (-1)^i * exp(-2 * i^2 * pi^2 *
                                               x)
    a
  }
  brp <- function(x, n = 1000) {
    a <- 1
    for (i in 1:n) a <- a - 2 * (-1)^(i - 1) * exp(-2 * i^2 *
                                                     x^2)
    a
  }
  global <- as.numeric(global & ncol(sresid) > 1)
  table <- NULL
  bbt <- as.list(1:(nvar + global))
  for (i in 1:nvar) {
    if (nvar != 1)
      usable <- which(varr[i, i, ] > 1e-12)
    else usable <- which(varr > 1e-12)
    w <- (n.risk[usable])^rho
    w <- w/sum(w)
    if (nvar != 1) {
      sci <- scaled[usable, i]
    }
    else sci <- scaled[usable]
    if (ties) {
      if (inherits(fit, "rsadd")) {
        sci <- as.vector(by(sci, time[usable], function(x) sum(x)/sqrt(length(x))))
        w <- as.vector(by(w, time[usable], sum))
      }
      else {
        w <- w * as.vector(table(time))[usable]
        w <- w/sum(w)
      }
    }
    sci <- sci * sqrt(w)
    timescale <- cumsum(w)
    bm <- cumsum(sci)
    bb <- bm - timescale * bm[length(bm)]
    if (test == "max")
      table <- rbind(table, c(max(abs(bb)), 1 - brp(max(abs(bb)))))
    else if (test == "cvm")
      table <- rbind(table, c(u2(bb), 1 - wc(u2(bb))))
    bbt[[i]] <- cbind(timescale, bb)
  }
  if (inherits(fit, "rsadd")) {
    beta <- fit$coef[1:(length(fit$coef) - length(fit$int) +  1)]
  }
  else beta <- fit$coef
  if (global) {
    qform <- function(matrix, vector) t(vector) %*% matrix %*%
      vector
    diagonal <- apply(varr, 3, diag)
    sumdiag <- apply(diagonal, 2, sum)
    usable <- which(sumdiag > 1e-12)
    score <- t(beta) %*% t(sresid[usable, ])
    varr <- varr[, , usable]
    qf <- apply(varr, 3, qform, vector = beta)
    w <- (n.risk[usable])^rho
    w <- w/sum(w)
    sci <- score/(qf)^0.5
    if (ties) {
      if (inherits(fit, "rsadd")) {
        sci <- as.vector(by(t(sci), time[usable], function(x) sum(x)/sqrt(length(x))))
        w <- as.vector(by(w, time[usable], sum))
      }
      else {
        w <- w * as.vector(table(time))
        w <- w/sum(w)
      }
    }
    sci <- sci * sqrt(w)
    timescale <- cumsum(w)
    bm <- cumsum(sci)
    bb <- bm - timescale * bm[length(bm)]
    if (test == "max")
      table <- rbind(table, c(max(abs(bb)), 1 - brp(max(abs(bb)))))
    else if (test == "cvm")
      table <- rbind(table, c(u2(bb), 1 - wc(u2(bb))))
    bbt[[nvar + 1]] <- cbind(timescale, bb)
    varnames <- c(varnames, "GLOBAL")
  }
  dimnames(table) <- list(varnames, c(test, "p"))
  out <- list(table = table, bbt = bbt, rho = rho)
  class(out) <- "rs.br"
  out
}



#' Behaviour of Covariates in Time for Relative Survival Regression Models
#'
#' Calculates the scaled partial residuals of a relative survival model
#' (\code{rsadd}, \code{rsmul} or \code{rstrans})
#'
#'
#' @param fit the result of fitting an additive relative survival model, using
#' the \code{rsadd}, \code{rsmul} or \code{rstrans} function.
#'
#' In the case of multiplicative and transformation models the output is
#' identical to \code{cox.zph} function, except no test is performed.
#' @param sc partial residuals calculated by the \code{resid} function. This is
#' used to save time if several tests are to be calculated on these residuals
#' and can otherwise be omitted.
#' @param transform a character string specifying how the survival times should
#' be transformed. Possible values are \code{"km"}, \code{"rank"},
#' \code{"identity"} and \code{log}. The default is \code{"identity"}.
#' @param var.type a character string specifying the variance used to scale the
#' residuals. Possible values are \code{"each"}, which estimates the variance
#' for each residual separately, and \code{sum}(default), which assumes the
#' same variance for all the residuals.
#' @return an object of class \code{rs.zph}. This function would usually be
#' followed by a plot of the result. The plot gives an estimate of the
#' time-dependent coefficient \code{beta(t)}. If the proportional hazards
#' assumption is true, \code{beta(t)} will be a horizontal line.
#' @seealso \code{\link{rsadd}}, \code{rsmul}, \code{rstrans},
#' \code{\link{resid}}, \code{\link[survival:cox.zph]{survival::cox.zph}}.
#' @references Goodness of fit: Stare J.,Pohar Perme M., Henderson R. (2005)
#' "Goodness of fit of relative survival models." Statistics in Medicine,
#' \bold{24}: 3911--3925.
#'
#' Package. Pohar M., Stare J. (2006) "Relative survival analysis in R."
#' Computer Methods and Programs in Biomedicine, \bold{81}: 272--278
#'
#' Relative survival: Pohar, M., Stare, J. (2007) "Making relative survival
#' analysis relatively easy." Computers in biology and medicine, \bold{37}:
#' 1741--1749.
#' @keywords survival
#' @examples
#'
#' data(slopop)
#' data(rdata)
#' fit <- rsadd(Surv(time,cens)~sex,rmap=list(age=age*365.241),
#' 	ratetable=slopop,data=rdata,int=5)
#' rszph <- rs.zph(fit)
#' plot(rszph)
#'
rs.zph <- function (fit, sc, transform = "identity", var.type = "sum")
{
  if (inherits(fit, "rsadd")) {
    if (missing(sc))
      sc <- resid(fit, "schoenfeld")
    sresid <- sc$res
    varr <- sc$kvarr
    fvar <- solve(sc$kvarr1)
    sresid <- as.matrix(sresid)
  }
  else {
    coef <- fit$coef
    options(warn = -1)
    sc <- coxph.detail(fit)
    options(warn = 0)
    sresid <- as.matrix(resid(fit, "schoenfeld"))
    varr <- sc$imat
    fvar <- fit$var
  }
  data <- fit$data[order(fit$data$Y), ]
  time <- data$Y
  stat <- data$stat

  if (!inherits(fit, "rsadd")) {
    ties <- as.vector(table(time[stat==1]))
    if(is.null(dim(varr))) varr <- rep(varr/ties,ties)
    else{
      varr <- apply(varr,1:2,function(x)rep(x/ties,ties))
      varr <- aperm(varr,c(2,3,1))
    }
  }
  keep <- 1:(length(fit$coef) - length(fit$int) + 1)
  varnames <- names(fit$coef)[keep]
  nvar <- length(varnames)
  ndead <- length(sresid)/nvar
  if (inherits(fit, "rsadd"))
    times <- time[stat == 1]
  else times <- sc$time
  if (is.character(transform)) {
    tname <- transform
    ttimes <- switch(transform, identity = times, rank = rank(times),
                     log = log(times), km = {
                       fity <- Surv(time, stat)
                       temp <- survfit(fity~1)
                       t1 <- temp$surv[temp$n.event > 0]
                       t2 <- temp$n.event[temp$n.event > 0]
                       km <- rep(c(1, t1), c(t2, 0))
                       if (is.null(attr(sresid, "strata")))
                         1 - km
                       else (1 - km[sort.list(sort.list(times))])
                     }, stop("Unrecognized transform"))
  }
  else {
    tname <- deparse(substitute(transform))
    ttimes <- transform(times)
  }
  if (var.type == "each") {
    invV <- apply(varr, 3, function(x) try(solve(x), silent = TRUE))
    if (length(invV) == length(varr)){
      if(!is.numeric(invV)){
        usable <- rep(FALSE, dim(varr)[3])
        options(warn=-1)
        invV <- as.numeric(invV)
        usable[1:(min(which(is.na(invV)))-1)] <- TRUE
        invV <- invV[usable]
        sresid <- sresid[usable,,drop=FALSE]
        options(warn=0)
      }
      else usable <- rep(TRUE, dim(varr)[3])
    }
    else {
      usable <- unlist(lapply(invV, is.matrix))
      if (!any(usable))
        stop("All the matrices are singular")
      invV <- invV[usable]
      sresid <- sresid[usable, , drop = FALSE]
    }
    di1 <- dim(varr)[1]
    di3 <- sum(usable)
    u <- array(data = matrix(unlist(invV), ncol = di1), dim = c(di1,
                                                                di1, di3))
    uv <- cbind(matrix(u, ncol = di1, byrow = TRUE), as.vector(t(sresid)))
    uv <- array(as.vector(t(uv)), dim = c(di1 + 1, di1, di3))
    r2 <- t(apply(uv, 3, function(x) x[1:di1, ] %*% x[di1 +
                                                        1, ]))
    r2 <- matrix(r2, ncol = di1)
    whr2 <-  apply(r2<100,1,function(x)!any(x==FALSE))
    usable <- as.logical(usable*whr2)
    r2 <- r2[usable,,drop=FALSE]
    u <- u[,,usable]
    dimnames(r2) <- list(times[usable], varnames)
    temp <- list(x = ttimes[usable], y = r2 + outer(rep(1,
                                                        sum(usable)), fit$coef[keep]), var = u, call = call,
                 transform = tname)
  }
  else if (var.type == "sum") {
    xx <- ttimes - mean(ttimes)
    r2 <- t(fvar %*% t(sresid) * ndead)
    r2 <- as.matrix(r2)
    dimnames(r2) <- list(times, varnames)
    temp <- list(x = ttimes, y = r2 + outer(rep(1, ndead),
                                            fit$coef[keep]), var = fvar, transform = tname)
  }
  else stop("Unknown 'var.type'")
  class(temp) <- "rs.zph"
  temp
}



#' Graphical Inspection of Proportional Hazards Assumption in Relative Survival
#' Models
#'
#' Displays a graph of the scaled partial residuals, along with a smooth curve.
#'
#'
#' @param x result of the \code{rs.zph} function.
#' @param resid a logical value, if \code{TRUE} the residuals are included on
#' the plot, as well as the smooth fit.
#' @param df the degrees of freedom for the fitted natural spline, \code{df=2}
#' leads to a linear fit.
#' @param nsmo number of points used to plot the fitted spline.
#' @param var the set of variables for which plots are desired.  By default,
#' plots are produced in turn for each variable of a model.  Selection of a
#' single variable allows other features to be added to the plot, e.g., a
#' horizontal line at zero or a main title.
#' @param cex a numerical value giving the amount by which plotting text and
#' symbols should be scaled relative to the default.
#' @param add logical, if \code{TRUE} the plot is added to an existing plot
#' @param col a specification for the default plotting color.
#' @param lty the line type.
#' @param xlab x axis label.
#' @param ylab y axis label.
#' @param xscale units for x axis, default is 1, i.e. days.
#' @param ... Additional arguments passed to the \code{plot} function.
#' @seealso \code{\link{rs.zph}}, \code{\link[survival:plot.cox.zph]{survival::plot.cox.zph}}
#' @references Goodness of fit: Stare J.,Pohar Perme M., Henderson R. (2005)
#' "Goodness of fit of relative survival models." Statistics in Medicine,
#' \bold{24}: 3911-3925.
#'
#' Package: Pohar M., Stare J. (2006) "Relative survival analysis in R."
#' Computer Methods and Programs in Biomedicine, \bold{81}: 272-278.
#'
#' Relative survival: Pohar, M., Stare, J. (2007) "Making relative survival
#' analysis relatively easy." Computers in biology and medicine, \bold{37}:
#' 1741-1749, 2007.
#' @keywords survival
#' @examples
#'
#' data(slopop)
#' data(rdata)
#' fit <- rsadd(Surv(time,cens)~sex+as.factor(agegr),rmap=list(age=age*365.241),
#'              ratetable=slopop,data=rdata,int=5)
#' rszph <- rs.zph(fit)
#' plot(rszph)
#'
plot.rs.zph <- function (x,resid = TRUE, df = 4, nsmo = 40, var, cex = 1,  add = FALSE, col = 1,
                         lty = 1, xlab, ylab, xscale = 1, ...)
{
  #require(splines)
  xx <- x$x
  if(x$transform=="identity")xx <- xx/xscale
  yy <- x$y
  d <- nrow(yy)
  df <- max(df)
  nvar <- ncol(yy)
  pred.x <- seq(from = min(xx), to = max(xx), length = nsmo)
  temp <- c(pred.x, xx)
  lmat <- splines::ns(temp, df = df, intercept = TRUE)
  pmat <- lmat[1:nsmo, ]
  xmat <- lmat[-(1:nsmo), ]
  qmat <- qr(xmat)
  if (missing(ylab))
    ylab <- paste("Beta(t) for", dimnames(yy)[[2]])
  if (missing(xlab))
    xlab <- "Time"
  if (missing(var))
    var <- 1:nvar
  else {
    if (is.character(var))
      var <- match(var, dimnames(yy)[[2]])
    if (any(is.na(var)) || max(var) > nvar || min(var) <
        1)
      stop("Invalid variable requested")
  }
  if (x$transform == "log") {
    xx <- exp(xx)
    pred.x <- exp(pred.x)
  }
  else if (x$transform != "identity") {
    xtime <- as.numeric(dimnames(yy)[[1]])/xscale
    apr1 <- approx(xx, xtime, seq(min(xx), max(xx), length = 17)[2 *
                                                                   (1:8)])
    temp <- signif(apr1$y, 2)
    apr2 <- approx(xtime, xx, temp)
    xaxisval <- apr2$y
    xaxislab <- rep("", 8)
    for (i in 1:8) xaxislab[i] <- format(temp[i])
  }
  for (i in var) {
    y <- yy[, i]
    yhat <- pmat %*% qr.coef(qmat, y)
    yr <- range(yhat, y)
    if (!add) {
      if (x$transform == "identity")
        plot(range(xx), yr, type = "n", xlab = xlab, ylab = ylab[i],...)
      else if (x$transform == "log")
        plot(range(xx), yr, type = "n", xlab = xlab, ylab = ylab[i],log = "x", ...)
      else {
        plot(range(xx), yr, type = "n", xlab = xlab, ylab = ylab[i],axes = FALSE, ...)
        axis(1, xaxisval, xaxislab)
        axis(2)
        box()
      }
    }
    if (resid)
      points(xx, y, cex = cex, col = col)
    lines(pred.x, yhat, col = col, lty = lty)
  }
}

plot.rs.br <- function (x, var, ylim = c(-2, 2), xlab, ylab, ...)
{
  bbt <- x$bbt
  par(ask = TRUE)
  if (missing(var))
    var <- 1:nrow(x$table)
  ychange <- FALSE
  if (missing(ylab))
    ylab <- paste("Brownian bridge for", row.names(x$table))
  else {
    if (length(ylab) == 1 & nrow(x$table) > 1)
      ylab <- rep(ylab, nrow(x$table))
  }
  if (missing(xlab))
    xlab <- "Time"
  for (i in var) {
    timescale <- bbt[[i]][, 1]
    bb <- bbt[[i]][, 2]
    plot(c(0, timescale), c(0, bb), type = "l", ylim = ylim,
         xlab = xlab, ylab = ylab[i], ...)
    abline(h = 1.36, col = 2)
    abline(h = 1.63, col = 2)
    abline(h = -1.36, col = 2)
    abline(h = -1.63, col = 2)
  }
  par(ask = FALSE)
}


Kernmatch <- function (t, tv, b, tD, nt4)
{
  kmat <- NULL
  for (it in 1:(length(nt4) - 1)) {
    kmat1 <- (outer(t[(nt4[it] + 1):nt4[it + 1]], tv, "-")/b[it])
    kmat1 <- kmat1^(kmat1 >= 0)
    kmat <- rbind(kmat, pmax(1 - kmat1^2, 0) * (1.5/b[it]))
  }
  kmat
}

kernerleftch <- function (td, b, nt4)
{
  n <- length(td)
  ttemp <- td[td >= b[1]]
  ntemp <- length(ttemp)
  if (ntemp == n)
    nt4 <- c(0, nt4[-1])
  else {
    nfirst <- n - ntemp
    nt4 <- c(0, 1:nfirst, nt4[-1])
    b <- c(td[1:nfirst], b)
  }
  krn <- Kernmatch(td, td, b, max(td), nt4)
  krn
}




#' Inverse transforming of time in Relative Survival
#'
#' This function can be used when predicting in Relative Survival using the
#' transformed time regression model (using \code{rstrans} function). It
#' inverses the time from Y to T in relative survival using the given
#' ratetable. The times Y can be produced with the \code{rstrans} function, in
#' which case, this is the reverse function. This function does the
#' transformation for one person at a time.
#'
#' Works only with ratetables that are split by age, sex and year. Transforming
#' can be computationally intensive, use lower and/or upper to guess the
#' interval of the result and thus speed up the function.
#'
#' @param y time in Y.
#' @param age age of the individual.  Must be in days.
#' @param sex sex of the individual. Must be coded in the same way as in the
#' \code{ratetable}.
#' @param year date of diagnosis. Must be in a Date or POSIXt format.
#' @param scale numeric value to scale the results. If \code{ratetable} is in
#' units/day, \code{scale = 365.241} causes the output to be reported in years.
#' @param ratetable a table of event rates, such as \code{survexp.us}.
#' @param lower the lower bound of interval where the result is expected. This
#' argument is optional, but, if given, can shorten the time the function needs
#' to calculate the result.
#' @param upper the upper bound of interval where the result is expected. See
#' \code{lower}.
#' @return A list of values \item{T}{the original time} \item{Y}{the
#' transformed time}.
#' @seealso \code{\link{rstrans}}
#' @references Package: Pohar M., Stare J. (2006) "Relative survival analysis
#' in R."  Computer Methods and Programs in Biomedicine, \bold{81}: 272-278.
#'
#' Relative survival: Pohar, M., Stare, J. (2007) "Making relative survival
#' analysis relatively easy." Computers in biology and medicine, \bold{37}:
#' 1741-1749.
#' @keywords survival
#' @examples
#'
#' data(slopop)
#' invtime(y = 0.1, age = 23011, sex = 1, year = as.Date('1986-01-01'), ratetable = slopop)
#'
invtime <- function (y = 0.1, age = 23011, sex = "male", year = as.Date('1986-01-01'), scale = 1,
                     ratetable = relsurv::slopop, lower, upper)
{
  if (!is.numeric(age))
    stop("\"age\" must be numeric", call. = FALSE)
  if (!is.numeric(y))
    stop("\"y\" must be numeric", call. = FALSE)
  if (!is.numeric(scale))
    stop("\"scale\" must be numeric", call. = FALSE)
  temp <- data.frame(age = age, sex = I(sex), year = year)
  if (missing(lower)) {
    if (!missing(upper))
      stop("Argument \"lower\" is missing, with no default",
           call. = FALSE)
    nyears <- round((110 - age/365.241))
    tab <- data.frame(age = rep(age, nyears), sex = I(rep(sex,
                                                          nyears)), year = rep(year, nyears))
    vred <- 1 - survexp(c(0, 1:(nyears - 1)) * 365.241 ~ ratetable(age = age,
                                                                   sex = sex, year = year), ratetable = ratetable, data = tab,
                        cohort = FALSE)
    place <- sum(vred <= y)
    if (place == 0)
      lower <- 0
    else lower <- floor((place - 1) * 365.241 - place)
    upper <- ceiling(place * 365.241 + place)
  }
  else {
    if (missing(upper))
      stop("Argument \"upper\" is missing, with no default",
           call. = FALSE)
    if (!is.integer(lower))
      lower <- floor(lower)
    if (!is.integer(upper))
      upper <- ceiling(upper)
    if (upper <= lower)
      stop("'upper' must be higher than 'lower'", call. = FALSE)
  }
  lower <- max(0, lower)
  tab <- data.frame(age = rep(age, upper - lower + 1), sex = I(rep(sex,
                                                                   upper - lower + 1)), year = rep(year, upper - lower +
                                                                                                     1))
  vred <- 1 - survexp((lower:upper) ~ ratetable(age = age,
                                                sex = sex, year = year), ratetable = ratetable, data = tab,
                      cohort = FALSE)
  place <- sum(vred <= y)
  if (place == 0)
    warning(paste("The event happened on or before day",
                  lower), call. = FALSE)
  if (place == length(vred))
    warning(paste("The event happened on or after day", upper),
            call. = FALSE)
  t <- (place + lower - 1)/scale
  age <- round(age/365.241, 0.01)
  return(list(age, sex, year, Y = y, T = t))
}






#' Fit Andersen et al Multiplicative Regression Model for Relative Survival
#'
#' Fits the Andersen et al multiplicative regression model in relative
#' survival. An extension of the coxph function using relative survival.
#'
#' NOTE: The follow-up time must be specified in days. The \code{ratetable}
#' being used may have different variable names and formats than the user's
#' data set, this is dealt with by the \code{rmap} argument. For example, if
#' age is in years in the data set but in days in the \code{ratetable} object,
#' age=age*365.241 should be used. The calendar year can be in any date format
#' (Date and POSIXt are allowed), the date formats in the
#' \code{ratetable} and in the data may differ.
#'
#' @param formula a formula object, with the response as a \code{Surv} object
#' on the left of a \code{~} operator, and, if desired, terms separated by the
#' \code{+} operator on the right.
#'
#' NOTE: The follow-up time must be in days.
#' @param data a data.frame in which to interpret the variables named in the
#' \code{formula}.
#' @param ratetable a table of event rates, such as \code{slopop}.
#' @param int the number of follow-up years used for calculating survival(the
#' data are censored after this time-point). If missing, it is set the the
#' maximum observed follow-up time.
#' @param na.action a missing-data filter function, applied to the model.frame,
#' after any subset argument has been used.  Default is
#' \code{options()$na.action}.
#' @param init vector of initial values of the iteration.  Default initial
#' value is zero for all variables.
#' @param method the default method \code{mul} assumes hazard to be constant on
#' yearly intervals. Method \code{mul1} uses the ratetable to determine the
#' time points when hazard changes. The \code{mul1} method is therefore more
#' accurate, but at the same time can be more computationally intensive.
#' @param control a list of parameters for controlling the fitting process.
#' See the documentation for \code{coxph.control} for details.
#' @param rmap an optional list to be used if the variables are not organized
#' and named in the same way as in the \code{ratetable} object. See details
#' below.
#' @param ... Other arguments will be passed to \code{coxph.control}.
#' @return an object of class \code{coxph} with an additional item:
#' \item{basehaz}{Cumulative baseline hazard (population values are seen as
#' offset) at centered values of covariates.}
#' @seealso \code{\link{rsadd}}, \code{\link{rstrans}}.
#' @references Method: Andersen, P.K., Borch-Johnsen, K., Deckert, T., Green,
#' A., Hougaard, P., Keiding, N. and Kreiner, S. (1985) "A Cox regression model
#' for relative mortality and its application to diabetes mellitus survival
#' data.", Biometrics, \bold{41}: 921--932.
#'
#' Package. Pohar M., Stare J. (2006) "Relative survival analysis in R."
#' Computer Methods and Programs in Biomedicine, \bold{81}: 272--278
#'
#' Relative survival: Pohar, M., Stare, J. (2007) "Making relative survival
#' analysis relatively easy." Computers in biology and medicine, \bold{37}:
#' 1741--1749.
#' @keywords survival
#' @examples
#'
#' data(slopop)
#' data(rdata)
#' #fit a multiplicative model
#' #note that the variable year is given in days since 01.01.1960 and that
#' #age must be multiplied by 365.241 in order to be expressed in days.
#' fit <- rsmul(Surv(time,cens)~sex+as.factor(agegr),rmap=list(age=age*365.241),
#'             ratetable=slopop,data=rdata)
#'
#'
#' #check the goodness of fit
#' rs.br(fit)
#'
#'
rsmul <- function (formula = formula(data), data = parent.frame(), ratetable = relsurv::slopop,
                   int, na.action, init, method = "mul", control,rmap, ...)
{
  #require(survival)

  if (!missing(rmap)) {
    rmap <- substitute(rmap)
  }
  rform <- rformulate(formula,data, ratetable, na.action,rmap,int)


  U <- rform$data
  if (missing(int))
    int <- ceiling(max(rform$Y/365.241))
  if(length(int)!=1)int <- max(int)
  fk <- (attributes(rform$ratetable)$factor != 1)
  nfk <- length(fk)
  if (method == "mul") {
    U <- survsplit(U, cut = (1:int) * 365.241, end = "Y",
                   event = "stat", start = "start", episode = "epi")
    fk <- (attributes(rform$ratetable)$factor != 1)
    nfk <- length(fk)
    U[, 4:(nfk + 3)] <- U[, 4:(nfk + 3)] + 365.241 * (U$epi) %*%
      t(fk)
    nsk <- dim(U)[1]
    xx <- exp_prep(U[, 4:(nfk + 3),drop=FALSE], 365.241, rform$ratetable)
    lambda <- -log(xx)/365.241
  }
  else if (method == "mul1") {
    U$id <- 1:dim(U)[1]
    my.fun <- function(x, attcut, nfk, fk) {
      intr <- NULL
      for (i in 1:nfk) {
        if (fk[i]) {
          n1 <- max(findInterval(as.numeric(x[3 + i]) +
                                   as.numeric(x[1]), attcut[[i]]) + 1, 2)
          n2 <- findInterval(as.numeric(x[3 + i]) + as.numeric(x[2]),
                             attcut[[i]])
          if (n2 > n1 & length(attcut[[i]] > 1)) {
            if (n2 > length(attcut[[i]]))
              n2 <- length(attcut[[i]])
            intr <- c(intr, as.numeric(attcut[[i]][n1:n2]) -
                        as.numeric(x[3 + i]))
          }
        }
      }
      intr <- sort(unique(c(intr, as.numeric(x[2]))))
      intr
    }
    attcut <- attributes(rform$ratetable)$cutpoints
    intr <- apply(U[, 1:(3 + nfk)], 1, my.fun, attcut, nfk,
                  fk)
    dolg <- unlist(lapply(intr, length))
    newdata <- lapply(U, rep, dolg)
    stoptime <- unlist(intr)
    starttime <- c(-1, stoptime[-length(stoptime)])
    first <- newdata$id != c(-1, newdata$id[-length(newdata$id)])
    starttime[first] <- newdata$start[first]
    last <- newdata$id != c(newdata$id[-1], -1)
    event <- rep(0, length(newdata$id))
    event[last] <- newdata$stat[last]
    U <- do.call("data.frame", newdata)
    U$start <- starttime
    U$Y <- stoptime
    U$stat <- event
    U[, 4:(nfk + 3)] <- U[, 4:(nfk + 3)] + (U$start) %*%
      t(fk)
    nsk <- dim(U)[1]
    xx <- exp_prep(U[, 4:(nfk + 3),drop=FALSE], 1, rform$ratetable)
    lambda <- -log(xx)/1
  }
  else stop("'method' must be one of 'mul' or 'mul1'")
  U$lambda <- log(lambda)
  if (rform$m == 0)
    fit <- coxph(Surv(start, Y, stat) ~ 1 + offset(lambda),
                 data = U, init = init, control = control, x = TRUE,
                 ...)
  else {
    xmat <- as.matrix(U[, (3 + nfk + 1):(ncol(U) - 2)])
    fit <- coxph(Surv(start, Y, stat) ~ xmat + offset(lambda),
                 data = U, init = init, control = control, x = TRUE,
                 ...)
    names(fit[[1]]) <- names(U)[(3 + nfk + 1):(ncol(U) -
                                                 2)]
  }
  class(fit) <- c("rsmul",class(fit))
  fit$basehaz <- basehaz(fit)			#NEW 2.05
  fit$data <- rform$data
  fit$call <- match.call()
  fit$int <- int
  if (length(rform$na.action))
    fit$na.action <- rform$na.action
  fit
}



#' Fit Cox Proportional Hazards Model in Transformed Time
#'
#' The function transforms each person's time to his/her probability of dying
#' at that time according to the ratetable. It then fits the Cox proportional
#' hazards model with the transformed times as a response. It can also be used
#' for calculatin the transformed times (no covariates are needed in the
#' formula for that purpose).
#'
#' NOTE: The follow-up time must be specified in days. The \code{ratetable}
#' being used may have different variable names and formats than the user's
#' data set, this is dealt with by the \code{rmap} argument. For example, if
#' age is in years in the data set but in days in the \code{ratetable} object,
#' age=age*365.241 should be used. The calendar year can be in any date format
#' (Date and POSIXt are allowed), the date formats in the
#' \code{ratetable} and in the data may differ.  A side product of this
#' function are the transformed times - stored in teh \code{y} object of the
#' output. To get these times, covariates are of course irrelevant.
#'
#' @param formula a formula object, with the response as a \code{Surv} object
#' on the left of a \code{~} operator, and, if desired, terms separated by the
#' \code{+} operator on the right.
#'
#' NOTE: The follow-up time must be in days.
#' @param data a data.frame in which to interpret the variables named in the
#' \code{formula}.
#' @param ratetable a table of event rates, such as \code{slopop}.
#' @param int the number of follow-up years used for calculating survival(the
#' rest is censored). If missing, it is set the the maximum observed follow-up
#' time.
#' @param na.action a missing-data filter function, applied to the model.frame,
#' after any subset argument has been used.  Default is
#' \code{options()$na.action}.
#' @param init vector of initial values of the iteration.  Default initial
#' value is zero for all variables.
#' @param control a list of parameters for controlling the fitting process.
#' See the documentation for \code{coxph.control} for details.
#' @param rmap an optional list to be used if the variables are not organized
#' and named in the same way as in the \code{ratetable} object. See details
#' below.
#' @param ... other arguments will be passed to \code{coxph.control}.
#' @return an object of class \code{coxph}. See \code{coxph.object} and
#' \code{coxph.detail} for details.  \item{y}{ an object of class \code{Surv}
#' containing the transformed times (these times do not depend on covariates).
#' }
#' @seealso \code{\link{rsmul}}, \code{\link{invtime}}, \code{\link{rsadd}},
#' \code{\link[survival:survexp]{survival::survexp}}.
#' @references Method: Stare J., Henderson R., Pohar M. (2005) "An individual
#' measure for relative survival." Journal of the Royal Statistical Society:
#' Series C, \bold{54} 115--126.
#'
#' Package. Pohar M., Stare J. (2006) "Relative survival analysis in R."
#' Computer Methods and Programs in Biomedicine, \bold{81}: 272--278
#'
#' Relative survival: Pohar, M., Stare, J. (2007) "Making relative survival
#' analysis relatively easy." Computers in biology and medicine, \bold{37}:
#' 1741--1749.
#' @keywords survival
#' @examples
#'
#' data(slopop)
#' data(rdata)
#'
#' #fit a Cox model using the transformed times
#' #note that the variable year is given in days since 01.01.1960 and that
#' #age must be multiplied by 365.241 in order to be expressed in days.
#' fit <- rstrans(Surv(time,cens)~sex+as.factor(agegr),rmap=list(age=age*365.241,
#'         sex=sex,year=year),ratetable=slopop,data=rdata)
#'
#'
#' #check the goodness of fit
#' rs.br(fit)
#'
rstrans <- function (formula = formula(data), data = parent.frame(), ratetable = relsurv::slopop,
                     int, na.action, init, control,rmap, ...)
{
  if (!missing(rmap)) {
    rmap <- substitute(rmap)
  }
  rform <- rformulate(formula, data, ratetable, na.action, rmap, int)


  if (missing(int))
    int <- ceiling(max(rform$Y/365.241))
  fk <- (attributes(rform$ratetable)$factor != 1)
  nfk <- length(fk)
  if (rform$type == "counting") {
    start <- 1 - exp_prep(rform$R, rform$start, rform$ratetable)
  }
  else start <- rep(0, rform$n)
  stop <- 1 - exp_prep(rform$R, rform$Y, rform$ratetable)
  if(any(stop==0&rform$Y!=0))stop[stop==0&rform$Y!=0] <- .Machine$double.eps
  if(length(int)!=1)int <- max(int)
  data <- rform$data
  stat <- rform$status
  if (rform$m == 0) {
    if (rform$type == "counting")
      fit <- coxph(Surv(start, stop, stat) ~ 1,
                   init = init, control = control, x = TRUE, ...)
    else fit <- coxph(Surv(stop, stat) ~ 1,
                      init = init, control = control, x = TRUE, ...)
  }
  else {
    xmat <- as.matrix(data[, (4 + nfk):ncol(data)])
    fit <- coxph(Surv(start, stop, stat) ~ xmat,
                 init = init, control = control, x = TRUE, ...)
    names(fit[[1]]) <- names(rform$X)
  }
  fit$call <- match.call()
  if (length(rform$na.action))
    fit$na.action <- rform$na.action
  data$start <- start
  data$Y <- stop
  fit$data <- data
  fit$int <- int
  return(fit)
}


#' Reorganize Data into a Ratetable Object
#'
#' The function assists in reorganizing certain types of data into a ratetable
#' object.
#'
#' This function only applies for ratetables that are organized by age, sex and
#' year.
#'
#' @param men a matrix containing the yearly (conditional) probabilities of one
#' year survival for men. Rows represent age (increasing 1 year per
#' line,starting with 0), the columns represent cohort years (the limits are in
#' \code{yearlim}, the increase is in \code{int.length}.
#' @param women a matrix containing the yearly (conditional) probabilities of
#' one year survival for women.
#' @param yearlim the first and last cohort year given in the tables.
#' @param int.length the length of intervals in which cohort years are given.
#' @return An object of class \code{ratetable}.
#' @seealso \code{\link[survival:ratetable]{survival::ratetable}}.
#' @references Package. Pohar M., Stare J. (2006) "Relative survival analysis
#' in R." Computer Methods and Programs in Biomedicine, \bold{81}: 272--278
#'
#' Relative survival: Pohar, M., Stare, J. (2007) "Making relative survival
#' analysis relatively easy." Computers in biology and medicine, \bold{37}:
#' 1741--1749.
#' @keywords survival
#' @examples
#'
#' men <- cbind(exp(-365.241*exp(-14.5+.08*(0:100))),exp(-365*exp(-14.7+.085*(0:100))))
#' women <- cbind(exp(-365.241*exp(-15.5+.085*(0:100))),exp(-365*exp(-15.7+.09*(0:100))))
#' table <- transrate(men,women,yearlim=c(1980,1990),int.length=10)
#'
transrate <- function (men, women, yearlim, int.length = 1)
{
  if (any(dim(men) != dim(women)))
    stop("The men and women matrices must be of the same size. \n In case of missing values at the end carry the last value forward")
  if ((yearlim[2] - yearlim[1])/int.length + 1 != dim(men)[2])
    stop("'yearlim' cannot be divided into intervals of equal length")
  if (!is.matrix(men) | !is.matrix(women))
    stop("input tables must be of class matrix")
  dimi <- dim(men)
  temp <- array(c(men, women), dim = c(dimi, 2))
  temp <- -log(temp)/365.241
  temp <- aperm(temp, c(1, 3, 2))
  cp <- as.Date(apply(matrix(yearlim[1] + int.length * (0:(dimi[2] -
                                                             1)), ncol = 1), 1, function(x) {
                                                               paste(x, "-01-01", sep = "")
                                                             }))
  attributes(temp) <- list(dim = c(dimi[1], 2, dimi[2]), dimnames = list(age=as.character(0:(dimi[1] -
                                                                                               1)), sex=c("male", "female"), year=as.character(yearlim[1] + int.length *
                                                                                                                                                 (0:(dimi[2] - 1)))), dimid = c("age", "sex", "year"),
                           factor = c(0, 1, 0),type=c(2,1,3), cutpoints = list((0:(dimi[1] - 1)) *
                                                                                 (365.241), NULL, cp), class = "ratetable")
  attributes(temp)$summary <- function (R)
  {
    x <- c(format(round(min(R[, 1])/365.241, 1)), format(round(max(R[,
                                                                     1])/365.241, 1)), sum(R[, 2] == 1), sum(R[, 2] == 2))
    x2 <- as.character(as.Date(c(min(R[, 3]), max(R[, 3])), origin=as.Date('1970-01-01')))
    paste("  age ranges from", x[1], "to", x[2], "years\n", " male:",
          x[3], " female:", x[4], "\n", " date of entry from",
          x2[1], "to", x2[2], "\n")
  }
  temp
}



#' Reorganize Data obtained from Human Life-Table Database into a Ratetable
#' Object
#'
#' The function assists in reorganizing the .txt files obtained from Human
#' Life-Table Database (http://www.lifetable.de -> Data by Country) into a
#' ratetable object.
#'
#' This function works with any table organised in the format provided by the
#' Human Life-Table Database, but currently only works with TypeLT 1 (i.e. age
#' intervals of length 1). The age must always start with value 0, but can end
#' at different values (when that happens, the last value is carried forward).
#' The rates between the cutpoints are taken to be constant.
#'
#' @param file a vector of file names which the data are to be read from. Must
#' be in .tex format and in the same format as the files in Human Life-Table
#' Database.
#' @param cut.year a vector of cutpoints for years. Must be specified when the
#' year spans in the files are not consecutive.
#' @param race a vector of race names for the input files.
#' @return An object of class \code{ratetable}.
#' @seealso \code{\link[survival:ratetable]{survival::ratetable}}, \code{\link{transrate.hmd}},
#' \code{\link{joinrate}}, \code{\link{transrate}}.
#' @references Package. Pohar M., Stare J. (2006) "Relative survival analysis
#' in R." Computer Methods and Programs in Biomedicine, \bold{81}: 272--278
#'
#' Relative survival: Pohar, M., Stare, J. (2007) "Making relative survival
#' analysis relatively easy." Computers in biology and medicine, \bold{37}:
#' 1741--1749.
#' @keywords survival
#' @examples
#'
#' \dontrun{
#' finpop <- transrate.hld(c("FIN_1981-85.txt","FIN_1986-90.txt","FIN_1991-95.txt"))
#' }
#' \dontrun{
#' nzpop <- transrate.hld(c("NZL_1980-82_Non-maori.txt","NZL_1985-87_Non-maori.txt",
#' 	 "NZL_1980-82_Maori.txt","NZL_1985-87_Maori.txt"),
#' 	 cut.year=c(1980,1985),race=rep(c("nonmaori","maori"),each=2))
#' }
#'
transrate.hld <- function(file, cut.year,race){
  nfiles <- length(file)
  data <- NULL
  for(it in 1:nfiles){
    tdata <- read.table(file[it],sep=",",header=TRUE)
    if(!any(tdata$TypeLT==1)) stop("Currently only TypeLT 1 is implemented")
    names(tdata) <- gsub(".","",names(tdata),fixed=TRUE)
    tdata <- tdata[,c("Country","Year1","Year2","TypeLT","Sex","Age","AgeInt","qx")]
    tdata <- tdata[tdata$TypeLT==1,]		#NEW - prej sem gledala tudi AgeInt, izkaze se, da ni treba. pri q(x) bi bilo vseeno tudi, ce bi gledala TypeLT=3.
    tdata <- tdata[!is.na(tdata$AgeInt),]		#NEW - vrzem ven zadnji interval, ki gre v neskoncnost in vsi umrejo (inf hazard)
    if(!missing(race))tdata$race <- rep(race[it],nrow(tdata))
    data <- rbind(data,tdata)
  }
  if(length(unique(data$Country))>1)warning("The data belongs to different countries")
  data <- data[order(data$Year1,data$Age),]
  data$qx <- as.character(data$qx)
  options(warn = -1)
  data$qx[data$qx=="."] <- NA
  data$qx <- as.numeric(data$qx)
  options(warn = 0)
  if(missing(cut.year)){
    y1 <-  unique(data$Year1)
    y2 <-  unique(data$Year2)
    if(any(apply(cbind(y1[-1],y2[-length(y2)]),1,diff)!=-1))warning("Data is not given for all the cut.year between the minimum and the maximum, use argument 'cut.year'")
  }
  else
    y1 <- cut.year
  if(length(y1)!=length(unique(data$Year1)))stop("Length 'cut.year' must match the number of unique values of Year1")
  cp <- as.Date(apply(matrix(y1,ncol=1),1,function(x){paste(x,"-01-01",sep="")}))
  dn2 <- as.character(y1)
  amax <- max(data$Age)
  a.fun <- function(data,amax){
    mdata <- data[data$Sex==1,]
    wdata <- data[data$Sex==2,]
    men <-NULL
    women <- NULL
    k <- sum(mdata$Age==0)
    mind <- c(which(mdata$Age[-nrow(mdata)] != mdata$Age[-1]-1),nrow(mdata))
    wind <- c(which(wdata$Age[-nrow(wdata)] != wdata$Age[-1]-1),nrow(wdata))
    mst <- wst <- 1
    for(it in 1:k){
      qx <- mdata[mst:mind[it],]$qx
      lqx <- length(qx)
      if(lqx!=amax+1){
        nmiss <- amax + 1 - lqx
        qx <- c(qx,rep(qx[lqx],nmiss))
      }
      naqx <- max(which(!is.na(qx)))
      if(naqx!=amax+1) qx[(naqx+1):(amax+1)] <- qx[naqx]
      men <- cbind(men,qx)
      mst <- mind[it]+1
      qx <- wdata[wst:wind[it],]$qx
      lqx <- length(qx)
      if(lqx!=amax+1){
        nmiss <- amax + 1 - lqx
        qx <- c(qx,rep(qx[lqx],nmiss))
      }
      naqx <- max(which(!is.na(qx)))
      if(naqx!=amax+1) qx[(naqx+1):(amax+1)] <- qx[naqx]
      women <- cbind(women,qx)
      wst <- wind[it]+1
    }
    men<- -log(1-men)/365.241
    women<- -log(1-women)/365.241
    dims <- c(dim(men),2)
    array(c(men,women),dim=dims)
  }
  if(missing(race)){
    out <- a.fun(data,amax)
    dims <- dim(out)
    attributes(out)<-list(
      dim=dims,
      dimnames=list(as.character(0:amax),as.character(y1),c("male","female")),
      dimid=c("age","year","sex"),
      factor=c(0,0,1),type=c(2,3,1),
      cutpoints=list((0:amax)*(365.241),cp,NULL),
      class="ratetable"
    )


  }
  else{
    race.val <- unique(race)
    if(length(race)!=length(file))stop("Length of 'race' must match the number of files")
    for(it in 1:length(race.val)){
      if(it==1){
        out <- a.fun(data[data$race==race.val[it],],amax)
        dims <- dim(out)
        out <- array(out,dim=c(dims,1))
      }
      else{
        out1 <- array(a.fun(data[data$race==race.val[it],],amax),dim=c(dims,1))
        out <- array(c(out,out1),dim=c(dims,it))
      }
    }
    attributes(out)<-list(
      dim=c(dims,it),
      dimnames=list(age=as.character(0:amax),year=as.character(y1),sex=c("male","female"),race=race.val),
      dimid=c("age","year","sex","race"),
      factor=c(0,0,1,1),type=c(2,3,1,1),
      cutpoints=list((0:amax)*(365.241),cp,NULL,NULL),
      class="ratetable"
    )
  }
  attributes(out)$summary <- function (R)
  {
    x <- c(format(round(min(R[, 1])/365.241, 1)), format(round(max(R[,
                                                                     1])/365.241, 1)), sum(R[, 3] == 1), sum(R[, 3] == 2))
    x2 <- as.character(as.Date(c(min(R[, 2]), max(R[, 2])), origin=as.Date('1970-01-01')))
    paste("  age ranges from", x[1], "to", x[2], "years\n", " male:",
          x[3], " female:", x[4], "\n", " date of entry from",
          x2[1], "to", x2[2], "\n")
  }
  out
}



#' Reorganize Data obtained from Human Mortality Database into a Ratetable
#' Object
#'
#' The function assists in reorganizing the .txt files obtained from Human
#' Mortality Database (http://www.mortality.org) into a ratetable object.
#'
#' This function works automatically with tables organised in the format
#' provided by the Human Mortality Database. Download Life Tables for Males and
#' Females separately from the column named 1x1 (period life tables, organized
#' by date of death, yearly cutpoints for age as well as calendar year).
#'
#' If you wish to provide the data in the required format by yourself, note
#' that the only two columns needed are calendar year (Year) and probability of
#' death (qx). Death probabilities must be calculated up to age 110 (in yearly
#' intervals).
#'
#' @param male a .txt file, containing the data on males.
#' @param female a .txt file, containing the data on females.
#' @return An object of class \code{ratetable}.
#' @seealso \code{\link[survival:ratetable]{survival::ratetable}}, \code{\link{transrate.hld}},
#' \code{\link{joinrate}}, \code{\link{transrate}}.
#' @references Package. Pohar M., Stare J. (2006) "Relative survival analysis
#' in R." Computer Methods and Programs in Biomedicine, \bold{81}: 272--278
#'
#' Relative survival: Pohar, M., Stare, J. (2007) "Making relative survival
#' analysis relatively easy." Computers in biology and medicine, \bold{37}:
#' 1741--1749.
#' @keywords survival
#' @examples
#'
#' \dontrun{
#' auspop <- transrate.hmd("mltper_1x1.txt","fltper_1x1.txt")
#' }
#'
transrate.hmd <- function(male,female){
  nfiles <- 2
  men <- try(read.table(male,sep="",header=TRUE),silent=TRUE)
  if(inherits(men, "try-error")){ men <- read.table(male,sep="",header=TRUE,skip=1)}
  men <- men[,c("Year","Age","qx")]
  y1 <- sort(unique(men$Year))
  ndata <- nrow(men)/111
  if(round(ndata)!=ndata)stop("Each year must contain ages from 0 to 110")
  men <- matrix(men$qx, ncol=ndata)
  men <- matrix(as.numeric(men),ncol=ndata)
  women <- try(read.table(female,sep="",header=TRUE),silent=TRUE)
  if(inherits(women, "try-error")) {women <- read.table(female,sep="",header=TRUE,skip=1)}
  women <- women[,"qx"]
  if(length(women)!=length(men))stop("Number of rows in the table must be equal for both sexes")
  women <- matrix(women, ncol=ndata)
  women <- matrix(as.numeric(women),ncol=ndata)
  cp <- as.Date(apply(matrix(y1, ncol = 1), 1, function(x) {
    paste(x, '-01-01', sep = "")
  }))
  dn2 <- as.character(y1)
  tfun <- function(vec){
    ind <- which(vec == 1 | is.na(vec))
    if(length(ind)>0)vec[min(ind):length(vec)] <- 0.999
    vec
  }
  men <- apply(men,2,tfun)
  women <- apply(women,2,tfun)
  men<- -log(1-men)/365.241
  women<- -log(1-women)/365.241
  nr <- nrow(men)-1
  dims <- c(dim(men),2)
  out <- array(c(men,women),dim=dims)
  attributes(out)<-list(
    dim=dims,
    dimnames=list(age=as.character(0:nr),year=as.character(y1),sex=c("male","female")),
    dimid=c("age","year","sex"),
    factor=c(0,0,1),type=c(2,3,1),
    cutpoints=list((0:nr)*(365.241),cp,NULL),
    class="ratetable"
  )
  attributes(out)$summary <- function (R)
  {
    x <- c(format(round(min(R[, 1])/365.241, 1)), format(round(max(R[,
                                                                     1])/365.241, 1)), sum(R[, 3] == 1), sum(R[, 3] == 2))
    x2 <- as.character(as.Date(c(min(R[, 2]), max(R[, 2])), origin=as.Date('1970-01-01')))
    paste("  age ranges from", x[1], "to", x[2], "years\n", " male:",
          x[3], " female:", x[4], "\n", " date of entry from",
          x2[1], "to", x2[2], "\n")
  }
  out
}






#' Join ratetables
#'
#' The function joins two or more objects organized as \code{ratetable} by
#' adding a new dimension.
#'
#' This function joins two or more \code{ratetable} objects by adding a new
#' dimension. The cutpoints of all the rate tables are compared; if merge=FALSE
#' (default) only the common intervals are kept, otherwise if merge=TRUE all
#' intervals are added (and hazards are extrapolated).
#' If the intervals defined by the cutpoints are not of
#' the same length, a warning message is displayed.  Each rate table must have
#' the same dimensions.
#'
#' @param tables a list of ratetables. If names are given, they are included as
#' \code{dimnames}.
#' @param dim.name the name of the added dimension.
#' @param merge if FALSE (default) only the intersect of all years/ages is taken.
#' If TRUE all possible years/ages are included (NOTE: in this case hazards are extrapolated
#' from earlier or later years/ages). This option only works for ratetables with dimensions
#  in the following order: age, year, sex.
#' @return An object of class \code{ratetable}.
#' @seealso \code{\link[survival:ratetable]{survival::ratetable}}, \code{\link{transrate.hld}},
#' \code{\link{transrate.hmd}}, \code{\link{transrate}}.
#' @references Package: Pohar M., Stare J. (2006) "Relative survival analysis
#' in R." Computer Methods and Programs in Biomedicine, \bold{81}: 272-278.
#'
#' Relative survival: Pohar, M., Stare, J. (2007) "Making relative survival
#' analysis relatively easy." Computers in biology and medicine, \bold{37}:
#' 1741-1749.
#' @keywords survival
#' @examples
#'
#' #newpop <- joinrate(list(Arizona=survexp.az,Florida=survexp.fl,
#' #                   Minnesota=survexp.mn),dim.name="state")
#'
joinrate <- function(tables,dim.name="country", merge=FALSE){
  nfiles <- length(tables)
  if(is.null(names(tables))) names(tables) <- paste("D",1:nfiles,sep="")
  if(any(!unlist(lapply(tables,is.ratetable))))stop("Tables must be in ratetable format")
  if(length(attributes(tables[[1]])$dim)<3)stop("Currently implemented for ratetables with 3 or more dimensions")

  if(is.null(attr(tables[[1]],"dimid")))attr(tables[[1]],"dimid") <- names((attr(tables[[1]],"dimnames")))

  if(merge){
    # Find all possible years and ages:
    add.years <- sort(unique(unlist(lapply(1:length(tables),
                                           function(x) attributes(tables[[x]])$dimnames$year))))
    add.ages <- as.character(sort(as.integer(unique(unlist(lapply(1:length(tables),
                                                                  function(x) attributes(tables[[x]])$dimnames$age))))))

    # Add the extra years/ages:
    for(tb in 1:length(tables)){
      tables[[tb]] <- ratetable_add_vals(tables[[tb]], add.years, add.ages)
    }
  }

  mc_length <- length(attributes(tables[[1]])$dimid)

  for(it in 2:nfiles){
    if(is.null(attr(tables[[it]],"dimid")))attr(tables[[it]],"dimid") <- names((attr(tables[[it]],"dimnames")))

    if(length(attributes(tables[[it]])$dimid)!=mc_length)stop("Each ratetable must have the same number of dimensions")
    mc <- match(attributes(tables[[it]])$dimid,attributes(tables[[1]])$dimid,nomatch=0)
    if(any(mc)==0) stop("Each ratetable must have the same number of dimensions")
    if(any(mc!=1:mc_length)){
      atts <- attributes(tables[[it]])
      tables[[it]] <- aperm(tables[[it]],mc)
      atts$dimid <- atts$dimid[mc]
      atts$dimnames <- atts$dimnames[mc]
      atts$cutpoints <- atts$cutpoints[mc]
      atts$factor <- atts$factor[mc]
      atts$type <- atts$type[mc]
      atts$dim <- atts$dim[mc]
      attributes(tables[[it]]) <- atts
    }
  }

  # Check if lists are equal:
  list.eq <- function(l1,l2){
    n <- length(l1)
    rez <- rep(TRUE,n)
    for(it in 1:n){
      if(length(l1[[it]])!=length(l2[[it]]))rez[it] <- FALSE
      else if(any(l1[[it]]!=l2[[it]]))rez[it] <- FALSE
    }
    rez
  }


  equal <- rep(TRUE,mc_length)
  for(it in 2:nfiles){
    equal <- equal*list.eq(attributes(tables[[1]])$cutpoints,attributes(tables[[it]])$cutpoints)
  }


  kir <-  which(!equal)

  newat <- attributes(tables[[1]])
  # imena <- list(d1=NULL,d2=NULL,d3=NULL,d4=NULL)
  imena <- eval(parse(text=paste0('list(', paste0('d', 1:mc_length, '=NULL', collapse = ', '), ')')))

  for(jt in kir){
    listy <- NULL
    for(it in 1:nfiles){
      listy <- c(listy,attributes(tables[[it]])$cutpoints[[jt]])
    }

    imena[[jt]] <- names(table(listy)[table(listy) == nfiles])

    # if(inherits(attributes(tables[[it]])$cutpoints[[jt]], 'date')){
    #   imena[[jt]] <- as.date(as.integer(imena[[jt]]))
    # }

    if(!length(imena[[jt]]))stop(paste("There are no common cutpoints for dimension", attributes(tables[[1]])$dimid[jt]))
  }

  # Merging:
  for(it in 1:nfiles){
    keep <- lapply(dim(tables[[it]]),function(x)1:x)
    for(jt in kir){
      meci <- which(match(attributes(tables[[it]])$cutpoints[[jt]],imena[[jt]],nomatch=0)!=0)

      if(it==1){
        newat$dimnames[[jt]] <- attributes(tables[[it]])$dimnames[[jt]][meci]
        newat$dim[[jt]] <- length(imena[[jt]])
        newat$cutpoints[[jt]] <- attributes(tables[[it]])$cutpoints[[jt]][meci]
      }
      if(length(meci)>1){if(max(diff(meci)!=1))warning(paste("The cutpoints for ",attributes(tables[[1]])$dimid[jt] ," are not equally spaced",sep=""))}
      keep[[jt]] <- meci
    }
    if(mc_length==3){
      tables[[it]] <- tables[[it]][keep[[1]],keep[[2]],keep[[3]]]
    } else if(mc_length==4){
      tables[[it]] <- tables[[it]][keep[[1]],keep[[2]],keep[[3]],keep[[4]]]
    } else if(mc_length==5){
      tables[[it]] <- tables[[it]][keep[[1]],keep[[2]],keep[[3]],keep[[4]],keep[[5]]]
    } else if(mc_length==6){
      tables[[it]] <- tables[[it]][keep[[1]],keep[[2]],keep[[3]],keep[[4]],keep[[5]],keep[[6]]]
    }
  }

  # Finish final object:
  dims <- newat$dim
  out <- array(tables[[1]],dim=c(dims,1))
  for(it in 2:nfiles){
    out1 <- array(tables[[it]],dim=c(dims,1))
    out <- array(c(out,out1),dim=c(dims,it))
  }
  mc <- 1:(mc_length+1)
  if(any(newat$factor>1)){
    wh <- which(newat$factor>1)
    mc <- c(mc[-wh],wh)
    out <- aperm(out,mc)
  }
  newat$dim <- c(dims,nfiles)[mc]
  newat$dimid <- c(newat$dimid,dim.name)[mc]
  if(mc_length==3){
    newat$cutpoints <- list(newat$cutpoints[[1]],newat$cutpoints[[2]],newat$cutpoints[[3]],NULL)[mc]
  } else if(mc_length==4){
    newat$cutpoints <- list(newat$cutpoints[[1]],newat$cutpoints[[2]],newat$cutpoints[[3]],newat$cutpoints[[4]],NULL)[mc]
  } else if(mc_length==5){
    newat$cutpoints <- list(newat$cutpoints[[1]],newat$cutpoints[[2]],newat$cutpoints[[3]],newat$cutpoints[[4]],newat$cutpoints[[5]],NULL)[mc]
  } else if(mc_length==6){
    newat$cutpoints <- list(newat$cutpoints[[1]],newat$cutpoints[[2]],newat$cutpoints[[3]],newat$cutpoints[[4]],newat$cutpoints[[5]],newat$cutpoints[[6]],NULL)[mc]
  }
  newat$factor <- c(newat$factor,1)[mc]
  newat$type <- c(newat$type,1)[mc]
  if(mc_length==3){
    newat$dimnames <- list(newat$dimnames[[1]],newat$dimnames[[2]],newat$dimnames[[3]],names(tables))[mc]
  } else if(mc_length==4){
    newat$dimnames <- list(newat$dimnames[[1]],newat$dimnames[[2]],newat$dimnames[[3]],newat$dimnames[[4]],names(tables))[mc]
  } else if(mc_length==5){
    newat$dimnames <- list(newat$dimnames[[1]],newat$dimnames[[2]],newat$dimnames[[3]],newat$dimnames[[4]],newat$dimnames[[5]],names(tables))[mc]
  } else if(mc_length==6){
    newat$dimnames <- list(newat$dimnames[[1]],newat$dimnames[[2]],newat$dimnames[[3]],newat$dimnames[[4]],newat$dimnames[[5]],newat$dimnames[[6]],names(tables))[mc]
  }
  names(newat$dimnames) <- newat$dimid
  attributes(out) <- newat
  out
}
# joinrate <- function(tables,dim.name="country", merge=FALSE){
#   nfiles <- length(tables)
#   if(is.null(names(tables))) names(tables) <- paste("D",1:nfiles,sep="")
#   if(any(!unlist(lapply(tables,is.ratetable))))stop("Tables must be in ratetable format")
#   if(length(attributes(tables[[1]])$dim)!=3)stop("Currently implemented only for ratetables with 3 dimensions")
#
#   if(is.null(attr(tables[[1]],"dimid")))attr(tables[[1]],"dimid") <- names((attr(tables[[1]],"dimnames")))
#
#   if(merge){
#     # Find all possible years and ages:
#     add.years <- sort(unique(unlist(lapply(1:length(tables),
#                                            function(x) attributes(tables[[x]])$dimnames$year))))
#     add.ages <- as.character(sort(as.integer(unique(unlist(lapply(1:length(tables),
#                                                                   function(x) attributes(tables[[x]])$dimnames$age))))))
#
#     # Add the extra years/ages:
#     for(tb in 1:length(tables)){
#       tables[[tb]] <- ratetable_add_vals(tables[[tb]], add.years, add.ages)
#     }
#   }
#
#   for(it in 2:nfiles){
#     if(is.null(attr(tables[[it]],"dimid")))attr(tables[[it]],"dimid") <- names((attr(tables[[it]],"dimnames")))
#     if(length(attributes(tables[[it]])$dimid)!=3)stop("Each ratetable must have 3 dimensions: age, year and sex")
#     mc <- match(attributes(tables[[it]])$dimid,attributes(tables[[1]])$dimid,nomatch=0)
#     if(any(mc)==0) stop("Each ratetable must have 3 dimensions: age, year and sex")
#     if(any(mc!=1:3)){
#       atts <- attributes(tables[[it]])
#       tables[[it]] <- aperm(tables[[it]],mc)
#       atts$dimid <- atts$dimid[mc]
#       atts$dimnames <- atts$dimnames[mc]
#       atts$cutpoints <- atts$cutpoints[mc]
#       atts$factor <- atts$factor[mc]
#       atts$type <- atts$type[mc]
#       atts$dim <- atts$dim[mc]
#       attributes(tables[[it]]) <- atts
#     }
#   }
#
#   # Check if lists are equal:
#   list.eq <- function(l1,l2){
#     n <- length(l1)
#     rez <- rep(TRUE,n)
#     for(it in 1:n){
#       if(length(l1[[it]])!=length(l2[[it]]))rez[it] <- FALSE
#       else if(any(l1[[it]]!=l2[[it]]))rez[it] <- FALSE
#     }
#     rez
#   }
#
#
#   equal <- rep(TRUE,3)
#   for(it in 2:nfiles){
#     equal <- equal*list.eq(attributes(tables[[1]])$cutpoints,attributes(tables[[it]])$cutpoints)
#   }
#
#
#   kir <-  which(!equal)
#
#   newat <- attributes(tables[[1]])
#   imena <- list(d1=NULL,d2=NULL,d3=NULL)
#
#   for(jt in kir){
#     listy <- NULL
#     for(it in 1:nfiles){
#       listy <- c(listy,attributes(tables[[it]])$cutpoints[[jt]])
#     }
#
#     imena[[jt]] <- names(table(listy)[table(listy) == nfiles])
#
#     if(inherits(attributes(tables[[it]])$cutpoints[[jt]], 'date')){
#       imena[[jt]] <- as.date(as.integer(imena[[jt]]))
#     }
#
#     if(!length(imena[[jt]]))stop(paste("There are no common cutpoints for dimension", attributes(tables[[1]])$dimid[jt]))
#   }
#
#   # Merging:
#   for(it in 1:nfiles){
#     keep <- lapply(dim(tables[[it]]),function(x)1:x)
#     for(jt in kir){
#       meci <- which(match(attributes(tables[[it]])$cutpoints[[jt]],imena[[jt]],nomatch=0)!=0)
#
#       if(it==1){
#         newat$dimnames[[jt]] <- attributes(tables[[it]])$dimnames[[jt]][meci]
#         newat$dim[[jt]] <- length(imena[[jt]])
#         newat$cutpoints[[jt]] <- attributes(tables[[it]])$cutpoints[[jt]][meci]
#       }
#       if(length(meci)>1){if(max(diff(meci)!=1))warning(paste("The cutpoints for ",attributes(tables[[1]])$dimid[jt] ," are not equally spaced",sep=""))}
#       keep[[jt]] <- meci
#     }
#     tables[[it]] <- tables[[it]][keep[[1]],keep[[2]],keep[[3]]]
#   }
#
#   # Finish final object:
#   dims <- newat$dim
#   out <- array(tables[[1]],dim=c(dims,1))
#   for(it in 2:nfiles){
#     out1 <- array(tables[[it]],dim=c(dims,1))
#     out <- array(c(out,out1),dim=c(dims,it))
#   }
#   mc <- 1:4
#   if(any(newat$factor>1)){
#     wh <- which(newat$factor>1)
#     mc <- c(mc[-wh],wh)
#     out <- aperm(out,mc)
#   }
#   newat$dim <- c(dims,nfiles)[mc]
#   newat$dimid <- c(newat$dimid,dim.name)[mc]
#   newat$cutpoints <- list(newat$cutpoints[[1]],newat$cutpoints[[2]],newat$cutpoints[[3]],NULL)[mc]
#   newat$factor <- c(newat$factor,1)[mc]
#   newat$type <- c(newat$type,1)[mc]
#   newat$dimnames <- list(newat$dimnames[[1]],newat$dimnames[[2]],newat$dimnames[[3]],names(tables))[mc]
#   names(newat$dimnames) <- newat$dimid
#   attributes(out) <- newat
#   out
# }



mlfit <- function (b, p, x, offset, d, h, ds, y, maxiter, tol)
{
  for (nit in 1:maxiter) {
    b0 <- b
    fd <- matrix(0, p, 1)
    sd <- matrix(0, p, p)
    if (nit == 1) {
      ebx <- exp(x %*% b) * exp(offset)
      l0 <- sum(d * log(h + ebx) - ds - y * ebx)
    }
    for (it in 1:p) {
      fd[it, 1] <- sum((d/(h + ebx) - y) * x[, it] * ebx)
      for (jt in 1:p) sd[it, jt] = sum((d/(h + ebx) - d *
                                          ebx/(h + ebx)^2 - y) * x[, it] * x[, jt] * ebx)
    }
    b <- b - solve(sd) %*% fd
    ebx <- exp(x %*% b) * exp(offset)
    l <- sum(d * log(h + ebx) - ds - y * ebx)
    bd <- abs(b - b0)
    if (max(bd) < tol)
      break()
  }
  out <- list(b = b, sd = sd, nit = nit, loglik = c(l0, l))
  out
}

print.rs.br <- function (x, digits = max(options()$digits - 4, 3), ...)
{
  invisible(print(x$table, digits = digits))
  if (x$rho != 0)
    invisible(cat("Weighted Brownian bridge with rho=", x$rho,
                  "\n"))
}

print.rsadd <- function (x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall: ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "", "\n")
  if (length(coef(x))) {
    cat("Coefficients")
    cat(":\n")
    print.default(format(x$coefficients, digits = digits),
                  print.gap = 2, quote = FALSE)
  }
  else cat("No coefficients\n\n")
  if(x$method=="EM")
    cat("\n", "Expected number of disease specific deaths: ",format(round(sum(x$Nie),2))," = ",format(round(100*sum(x$Nie)/sum(x$data$stat),1)),"% \n" ,sep="")
  if(x$method=="EM"|x$method=="max.lik"){
    chi <- 2*max((x$loglik[2]-x$loglik[1]),0)
    if(x$method=="EM")df <- length(x$coef)
    else df <- length(x$coef)-length(x$int)+1
    if(df>0){
      p.val <- 1- pchisq(chi,df)
      if(x$method=="max.lik")cat("\n")
      cat("Likelihood ratio test=",format(round(chi,2)),", on ",df," df, p=",format(p.val),"\n",sep="")
    }
    else cat("\n")
  }
  cat("n=",nrow(x$data),sep="")
  if(length(x$na.action))cat("  (",length(x$na.action)," observations deleted due to missing)",sep="")
  cat("\n")
  if (length(x$warnme))
    cat("\n", x$warnme, "\n\n")
  else cat("\n")
  invisible(x)
}

summary.rsadd <- function (object, correlation = FALSE, symbolic.cor = FALSE,
                           ...)
{
  if (inherits(object, "glm")) {
    p <- object$rank
    if (p > 0) {
      p1 <- 1:p
      Qr <- object$qr
      aliased <- is.na(coef(object))
      coef.p <- object$coefficients[Qr$pivot[p1]]
      covmat <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
      dimnames(covmat) <- list(names(coef.p), names(coef.p))
      var.cf <- diag(covmat)
      s.err <- sqrt(var.cf)
      tvalue <- coef.p/s.err
      dn <- c("Estimate", "Std. Error")
      pvalue <- 2 * pnorm(-abs(tvalue))
      coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
      dimnames(coef.table) <- list(names(coef.p), c(dn,
                                                    "z value", "Pr(>|z|)"))
      df.f <- NCOL(Qr$qr)
    }
    else {
      coef.table <- matrix(, 0, 4)
      dimnames(coef.table) <- list(NULL, c("Estimate",
                                           "Std. Error", "t value", "Pr(>|t|)"))
      covmat.unscaled <- covmat <- matrix(, 0, 0)
      aliased <- is.na(coef(object))
      df.f <- length(aliased)
    }
    ans <- c(object[c("call", "terms", "family", "iter",
                      "warnme")], list(coefficients = coef.table, var = covmat,
                                       aliased = aliased))
    if (correlation && p > 0) {
      dd <- s.err
      ans$correlation <- covmat/outer(dd, dd)
      ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.rsadd"
  }
  else if (inherits(object, "rsadd")) {
    aliased <- is.na(coef(object))
    coef.p <- object$coef
    var.cf <- diag(object$var)
    s.err <- sqrt(var.cf)
    tvalue <- coef.p/s.err
    dn <- c("Estimate", "Std. Error")
    pvalue <- 2 * pnorm(-abs(tvalue))
    coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
    dimnames(coef.table) <- list(names(coef.p), c(dn, "z value",
                                                  "Pr(>|z|)"))
    ans <- c(object[c("call", "terms", "iter", "var")], list(coefficients = coef.table,
                                                             aliased = aliased))
    if (correlation && sum(aliased) != length(aliased)) {
      dd <- s.err
      ans$correlation <- object$var/outer(dd, dd)
      ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.rsadd"
  }
  else ans <- object
  return(ans)
}

print.summary.rsadd <- function (x, digits = max(3, getOption("digits") - 3), symbolic.cor = x$symbolic.cor,
                                 signif.stars = getOption("show.signif.stars"), ...)
{
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  if (length(x$aliased) == 0) {
    cat("\nNo Coefficients\n")
  }
  else {
    cat("\nCoefficients:\n")
    coefs <- x$coefficients
    if (!is.null(aliased <- x$aliased) && any(aliased)) {
      cn <- names(aliased)
      coefs <- matrix(NA, length(aliased), 4, dimnames = list(cn,
                                                              colnames(coefs)))
      coefs[!aliased, ] <- x$coefficients
    }
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                 na.print = "NA", ...)
  }
  if (length(x$warnme))
    cat("\n", x$warnme, "\n")
  correl <- x$correlation
  if (!is.null(correl)) {
    p <- NCOL(correl)
    if (p > 1) {
      cat("\nCorrelation of Coefficients:\n")
      if (is.logical(symbolic.cor) && symbolic.cor) {
        print(symnum(correl, abbr.colnames = NULL))
      }
      else {
        correl <- format(round(correl, 2), nsmall = 2,
                         digits = digits)
        correl[!lower.tri(correl)] <- ""
        print(correl[-1, -p, drop = FALSE], quote = FALSE)
      }
    }
  }
  cat("\n")
  invisible(x)
}



#' Excess hazard function smoothing
#'
#' An Epanechnikov kernel function based smoother for smoothing the baseline
#' excess hazard calculated by the \code{rsadd} function with the \code{EM}
#' method.
#'
#' The function performs Epanechnikov kernel smoothing. The follow up time is
#' divided (according to percentiles of event times) into several intervals
#' (number of intervals defined by \code{n.bwin}) in which the width is
#' calculated as a factor of the maximum span between event times.  Boundary
#' effects are also taken into account on both sides.
#'
#' @param fit Fit from the additive relative survival model using the \code{EM}
#' method.
#' @param bwin The relative width of the smoothing window (default is 1).
#' @param times The times at which the smoother is to be evaluated. If missing,
#' it is evaluated at all event times.
#' @param n.bwin Number of times that the window width may change.
#' @param left If \code{FALSE} (default) smoothing is performed symmetrically,
#' if \code{TRUE} only leftside neighbours are considered.
#' @return A list with two components: \item{lambda}{the smoothed excess
#' baseline hazard function} \item{times}{the times at which the smoothed
#' excess baseline hazard is evaluated.}
#' @seealso \code{\link{rsadd}},
#' @references Package. Pohar M., Stare J. (2006) "Relative survival analysis
#' in R." Computer Methods and Programs in Biomedicine, \bold{81}: 272--278
#'
#' Relative survival: Pohar, M., Stare, J. (2007) "Making relative survival
#' analysis relatively easy." Computers in biology and medicine, \bold{37}:
#' 1741--1749.
#'
#' EM algorithm: Pohar Perme M., Henderson R., Stare, J. (2009) "An approach to
#' estimation in relative survival regression." Biostatistics, \bold{10}:
#' 136--146.
#' @keywords survival
#' @examples
#'
#' data(slopop)
#' data(rdata)
#' #fit an additive model with the EM method
#' fit <- rsadd(Surv(time,cens)~sex+age,rmap=list(age=age*365.241),
#' 		ratetable=slopop,data=rdata,int=5,method="EM")
#' sm <- epa(fit)
#' plot(sm$times,sm$lambda)
#'
epa <- function(fit,bwin,times,n.bwin=16,left=FALSE){
  #bwin ... width of the window, relative to the default (1)
  #fit ... EM fit
  #times... times at which the smoothed plot is calculated
  #n.bwin ... number of different windows
  #left ... only predictable smoothing
  utd <- fit$times
  if(missing(times))times <- seq(1,max(utd),length=100)
  if(max(times)>max(utd)){
    warning("Cannot extrapolate beyond max event time")
    times <- pmax(times,max(utd))
  }
  nutd <- length(utd)
  nt4 <- c(1,ceiling(nutd*(1:n.bwin)/n.bwin))
  if(missing(bwin))bwin <- rep(length(fit$times)/100,n.bwin)
  else bwin <- rep(bwin*length(fit$times)/100,n.bwin)
  for(it in 1:n.bwin){
    bwin[it] <- bwin[it]*max(diff(utd[nt4[it]:nt4[it+1]]))
  }
  while(utd[nt4[2]]<bwin[1]){		# ce je bwin velik, skrajsamo nt4
    nt4 <- nt4[-2]
    if(length(nt4)==1)break
  }
  #the smoothing matrix
  if(left) krn <- kernerleftch(utd,bwin,nt4)
  else krn <- kern(times,utd,bwin,nt4)
  lams <- pmax(krn%*%fit$lam0.ns,0)
  list(lambda=lams,times=times)				#	, weights=c(fit$times[1],diff(fit$times)))
}

Kern <- function (t, tv, b, tD, nt4)
{
  Rb <- max(tv)					#Right border
  kmat <- NULL
  tvs <- tv
  tv <- tv[-1]
  kt <- function(q,t)12*(t+1)/(1+q)^4*( (1-2*q)*t + (3*q^2-2*q+1)/2  )
  totcajti <- NULL
  for (it in 1:(length(nt4) - 1)) {
    cajti <- t[t>tvs[nt4[it]] & t<=tvs[nt4[it + 1]]]
    if(length(cajti)){
      q <- min( cajti/b[it],1,(Rb-cajti)/b[it])
      if(q<1 & length(cajti)>1){
        jc <- 1
        while(jc <=length(cajti)){
          qd <- pmin( cajti[jc:length(cajti)]/b[it],1,(Rb-cajti[jc:length(cajti)])/b[it])
          q <- qd[1]
          if(q==1){
            casi <- cajti[jc:length(cajti)][qd==1]
            q <- 1
            jc <- sum(qd==1)+jc
          }
          else{
            casi <- cajti[jc]
            jc <- jc+1
          }
          kmat1 <- outer(casi, tv, "-")/b[it]		#z - to je ok
          if(q<1){
            if(casi>b[it]) kmt1 <- -kmat1
            vr <- kt(q,kmat1)*(kmat1>=-1 & kmat1 <= q)
          }
          else vr <-  pmax((1 - kmat1^2) * .75,0)
          kmat <- rbind(kmat, vr/b[it])
          totcajti <- c(totcajti,casi)
        }
      }
      else{
        kmat1 <- outer(cajti, tv, "-")/b[it]		#z - to je ok
        q <- min( cajti/b[it],1)
        if(q<1)vr <- kt(q,kmat1)*(kmat1>=-1 & kmat1 <= q)
        else vr <-  pmax((1 - kmat1^2) * .75,0)
        kmat <- rbind(kmat, vr/b[it])
        totcajti <- c(totcajti,cajti)
      }#else
    }#if

  }#for
  kmat
}

kern <- function (times,td, b, nt4)
{
  n <- length(td)
  ttemp <- td[td >= b[1]]
  ntemp <- length(ttemp)
  if (ntemp == n)
    nt4 <- c(0, nt4[-1])
  td <- c(0,td)
  nt4 <- c(1,nt4+1)
  b <- c(b[1],b)
  krn <- Kern(times, td, b, max(td), nt4)
  krn
}

exp_prep <- function (x, y,ratetable,status,times,fast=FALSE,ys,prec,cmp=F,netweiDM=FALSE) {			#function that prepares the data for C function call

  #x= matrix of demographic covariates - each individual has one line
  #y= follow-up time for each individual (same length as nrow(x)!)
  #ratetable= rate table used for calculation
  #status= status for each individual (same length as nrow(x)!), not needed if we only need Spi, status needed for rs.surv
  #times= times at which we wish to evaluate the quantities, not needed if we only need Spi, times needed for rs.surv
  #fast=for mpp method only
  #netweiDM=should new netwei script be used

  x <- as.matrix(x)
  if (ncol(x) != length(dim(ratetable)))
    stop("x matrix does not match the rate table")
  atts <- attributes(ratetable)

  cuts <- atts$cutpoints

  if (is.null(atts$type)) {
    rfac <- atts$factor
    us.special <- (rfac > 1)
  }
  else {
    rfac <- 1 * (atts$type == 1)
    us.special <- (atts$type == 4)
  }
  if (length(rfac) != ncol(x))
    stop("Wrong length for rfac")


  if (any(us.special)) {
    if (sum(us.special) > 1)
      stop("Two columns marked for special handling as a US rate table")
    cols <- match(c("age", "year"), atts$dimid)
    if (any(is.na(cols)))
      stop("Ratetable does not have expected shape")
    if (exists("as.Date")) {
      bdate <- as.Date("1960/1/1") + (x[, cols[2]] - x[,
                                                       cols[1]])
      byear <- format(bdate, "%Y")
      offset <- as.numeric(bdate - as.Date(paste(byear,
                                                 "01/01", sep = "/")))
    } else stop("Can't find an appropriate date class\n")
    # Prej je bilo tole:
    # else if (exists("date.mdy")) {
    #   bdate <- as.date(x[, cols[2]] - x[, cols[1]])
    #   byear <- date.mdy(bdate)$year
    #   offset <- bdate - mdy.date(1, 1, byear)
    # }

    x[, cols[2]] <- x[, cols[2]] - offset
    if (any(rfac > 1)) {
      temp <- which(us.special)
      nyear <- length(cuts[[temp]])
      nint <- rfac[temp]
      cuts[[temp]] <- round(approx(nint * (1:nyear), cuts[[temp]],
                                   nint:(nint * nyear))$y - 1e-04)
    }
  }

  if(!missing(status)){		#the function was called from rs.surv
    if(length(status)!=nrow(x))    stop("Wrong length for status")

    if(missing(times))    times <- sort(unique(y))

    if (any(times < 0))
      stop("Negative time point requested")
    ntime <- length(times)
    if(missing(ys)) ys <- rep(0,length(y))
    #    times2 <- times
    #    times2[1] <- preci

    # It may be an integer...check and make it numeric. Otherwise problems in C functions:
    if(is.integer(x[1,1])){
      x <- apply(x, 2, as.numeric)
    }

    if(cmp)   temp <- .Call("cmpfast",  as.integer(rfac), 		#fast=pohar-perme or ederer2 - data from pop. tables only while under follow-up
                            as.integer(atts$dim), as.double(unlist(cuts)), ratetable,
                            x, y, ys,as.integer(status), times,PACKAGE="relsurv")
    else if(fast&!missing(prec))    temp <- .Call("netfastpinter2",  as.integer(rfac), 		#fast=pohar-perme or ederer2 - data from pop. tables only while under follow-up
                                                  as.integer(atts$dim), as.double(unlist(cuts)), ratetable,
                                                  x, y, ys,as.integer(status), times,prec,PACKAGE="relsurv")
    else if(fast&missing(prec))    temp <- .Call("netfastpinter",  as.integer(rfac), 		#fast=pohar-perme or ederer2 - data from pop. tables only while under follow-up
                                                 as.integer(atts$dim), as.double(unlist(cuts)), ratetable,
                                                 x, y, ys,as.integer(status), times,PACKAGE="relsurv")
    else if(netweiDM==TRUE)    temp <- .Call("netweiDM",  as.integer(rfac),
                                             as.integer(atts$dim), as.double(unlist(cuts)), ratetable,
                                             x, y, ys,as.integer(status), times,PACKAGE="relsurv")
    else    temp <- .Call("netwei",  as.integer(rfac),
                          as.integer(atts$dim), as.double(unlist(cuts)), ratetable,
                          x, y, as.integer(status), times,PACKAGE="relsurv")
  }
  else{				#only expected survival at time y is needed for each individual
    if(length(y)==1)y <- rep(y,nrow(x))
    if(length(y)!=nrow(x)) stop("Wrong length for status")
    temp <- .Call("expc",  as.integer(rfac),
                  as.integer(atts$dim), as.double(unlist(cuts)), ratetable,
                  x, y,PACKAGE="relsurv")
    temp  <- temp$surv
  }
  temp
}



#' Compute a Relative Survival Curve
#'
#' Computes an estimate of the relative survival curve using the Ederer I,
#' Ederer II method, Pohar-Perme method or the Hakulinen method
#'
#' NOTE: The follow-up time must be specified in days. The \code{ratetable}
#' being used may have different variable names and formats than the user's
#' data set, this is dealt with by the \code{rmap} argument. For example, if
#' age is in years in the data set but in days in the \code{ratetable} object,
#' age=age*365.241 should be used. The calendar year can be in any date format
#' (Date and POSIXt are allowed), the date formats in the
#' \code{ratetable} and in the data may differ.
#'
#' The potential censoring times needed for the calculation of the expected
#' survival by the Hakulinen method are calculated automatically. The times of
#' censoring are left as they are, the times of events are replaced with
#' \code{fin.date - year}.
#'
#' The calculation of the Pohar-Perme estimate is more time consuming since
#' more data are needed from the population tables.  The old version of the
#' function, now named \code{rs.survo} can be used as a faster version for the
#' Hakulinen and Ederer II estimate.
#'
#' Numerical integration is required for Pohar-Perme estimate. The integration
#' precision is set with argument \code{precision}, which defaults to daily
#' intervals, a default that should give enough precision for any practical
#' purpose.
#'
#' Note that even though the estimate is always calculated using numerical
#' integration, only the values at event and censoring times are reported.
#' Hence, the function \code{plot} draws a step function in between and the
#' function \code{summary} reports the value at the last event or censoring
#' time before the specified time. If the output of the estimated values at
#' other points is required, this should be specified with argument
#' \code{add.times}.
#'
#' Standardized net survival can be also calculated:
#'
#' \eqn{SNS (t) = \sum_{j} w_j {NS}_j (t)}
#'
#' where \eqn{NS_j} is the net survival in the j-th group, \eqn{w_j}
#' is the weight for the j-th group and \eqn{SNS} is the standardized net survival.
#' \eqn{SNS} can be calculated by using the \code{weight.table} and \code{weight.names} arguments. The function
#' also returns the corresponding standard error and confidence interval.
#'
#' @param formula a formula object, with the response as a \code{Surv} object
#' on the left of a \code{~} operator, and, if desired, terms separated by the
#' \code{+} operator on the right. If no strata are used, \code{~1} should be
#' specified. Allows for left-truncation, i.e. Surv(start,stop,event) formula.
#'
#' NOTE: The follow-up time must be in days.
#' @param data a data.frame in which to interpret the variables named in the
#' \code{formula}.
#' @param ratetable a table of event rates, organized as a \code{ratetable}
#' object, such as \code{slopop}.
#' @param na.action a missing-data filter function, applied to the model.frame,
#' after any subset argument has been used.  Default is
#' \code{options()$na.action}.
#' @param fin.date the date of the study ending, used for calculating the
#' potential follow-up times in the Hakulinen method. If missing, it is
#' calculated as \code{max(year+time)}.
#' @param method the method for calculating the relative survival. The options
#' are \code{pohar-perme}(default), \code{ederer1}, \code{ederer2} and
#' \code{hakulinen}.
#' @param conf.type one of \code{plain}, \code{log} (the default), or
#' \code{log-log}. The first option causes the standard intervals curve +- k
#' *se(curve), where k is determined from conf.int. The log option calculates
#' intervals based on the cumulative hazard or log(survival). The last option
#' bases intervals on the log hazard or log(-log(survival)).
#' @param conf.int the level for a two-sided confidence interval on the
#' survival curve(s). Default is 0.95.
#' @param type defines how survival estimates are to be calculated given the
#' hazards.  The default (\code{kaplan-meier}) calculates the product integral,
#' whereas the option \code{fleming-harrington} exponentiates the negative
#' cumulative hazard. Analogous to the usage in \code{survfit}.
#' @param add.times specific times at which the curve should be evaluated.
#' @param precision Precision for numerical integration. Default is 1, which
#' means that daily intervals are taken, the value may be decreased to get a
#' higher precision or increased to achieve a faster calculation. The
#' calculation intervals always include at least all times of event and
#' censoring as border points.
#' @param rmap an optional list to be used if the variables are not organized
#' and named in the same way as in the \code{ratetable} object. See details
#' below.
#' @param weight.table Default is NULL. If defined it is used for calculating
#' standardized net survival. Supply a data.frame
#' that contains the weights and group names for calculating the standardized net survival.
#' The group column also has to be included in the data argument.
#' For the theory see the details below.
#' @param weight.names Default is NULL. If defined it is used for calculating
#' standardized net survival. Supply a character vector of length two with the
#' names of the group and weight columns in \code{weight.table}.
#' @return a \code{survfit} object; see the help on \code{survfit.object} for
#' details.  The \code{survfit} methods are used for \code{print},
#' \code{summary}, \code{plot}, \code{lines}, and \code{points}.
#' @seealso \code{survfit}, \code{survexp}
#' @references Package: Pohar Perme, M., Pavlic, K. (2018) "Nonparametric
#' Relative Survival Analysis with the R Package relsurv". Journal of
#' Statistical Software. 87(8), 1-27, doi: "10.18637/jss.v087.i08" Theory:
#' Pohar Perme, M., Esteve, J., Rachet, B. (2016) "Analysing Population-Based
#' Cancer Survival - Settling the Controversies." BMC Cancer, 16 (933), 1-8.
#' doi:10.1186/s12885-016-2967-9. Theory: Pohar Perme, M., Stare, J., Esteve,
#' J. (2012) "On Estimation in Relative Survival", Biometrics, 68(1), 113-120.
#' doi:10.1111/j.1541-0420.2011.01640.x.
#' @keywords survival
#' @examples
#'
#' data(slopop)
#' data(rdata)
#' #calculate the relative survival curve
#' #note that the variable year must be given in a date format and that
#' #age must be multiplied by 365.241 in order to be expressed in days.
#' rs.surv(Surv(time,cens)~sex,rmap=list(age=age*365.241), ratetable=slopop,data=rdata)
#'
#' # SNS:
#' wei <- data.frame(agegr=c('<54', '54-61', '62-70', '71-95'),
#'                   weight=c(0.2, 0.2, 0.3, 0.3))
#' rs.surv(Surv(time,cens)~1, rmap=list(age=age*365.241), ratetable=slopop, data=rdata,
#'    weight.table=wei, weight.names = c('agegr', 'weight'))
#'
rs.surv <- function(formula = formula(data), data = parent.frame(),ratetable = relsurv::slopop,
                    na.action, fin.date, method = "pohar-perme", conf.type = "log",
                    conf.int = 0.95,type="kaplan-meier",add.times,precision=1,rmap,
                    weight.table=NULL, weight.names=NULL){

  # SNS part:
  if(!is.null(weight.table) & (!is.numeric(weight.table))){
    call <- match.call()

    # MAKE ALL CHECKS:
    if(!(inherits(weight.table, 'data.frame')
       # & all(colnames(weight.table) %in% c('weight', 'group'))
       # & ncol(weight.table) == 2
       )){
      stop("In the argument weight.table you should supply a data.frame.")
    }

    if(!(inherits(weight.names, 'character'))){
      stop("In the argument weight.names you should supply a character vector of length 2 containing the names of the group and weight columns in weight.table.")
    }

    if(!all(weight.names %in% colnames(weight.table))){
      stop('The table in weight.table does not contain the columns given in weight.names.')
    }

    if(!(weight.names[1] %in% colnames(data))){
      stop('The data does not contain the grouping column (first value) given in weight.names.')
    }

    group_val <- weight.table[,weight.names[1]]
    weight_val <- weight.table[,weight.names[2]]

    if(!all(group_val %in% unique(data[,weight.names[1]]))){
      stop("Not all weight categories from weight.table are present in the data.")
    }

    if(!all(unique(data[,weight.names[1]]) %in% group_val)){
      stop("Not all weight categories from the data are present in weight.table.")
    }

    if(abs(sum(weight_val) - 1) > 1e-16){
      # weight.table$weight <- weight.table$weight/sum(weight.table$weight)
      stop('The sum of the weights in the weight.table argument is not equal to 1.')
    }

    # Find groups for SNS:
    groups_u <- unique(group_val)
    groups_l <- length(groups_u)

    # Fix rmap:
    if (!missing(rmap)) {
      rmap <- substitute(rmap)
    }

    # Get full model:
    mod0 <- rs.surv(formula, data, ratetable,
                    na.action, fin.date, method, conf.type,
                    conf.int, type, add.times, precision=1, rmap,
                    weight.table=-999) # add weights arg so that it is detected later for rmap

    # Prepare object for output:
    out <- mod0

    # Save separate models:
    mod_l <- list()

    # Save maximum follow-up times among groups:
    if('strata' %in% names(mod0)){
      max.times <- data.frame(matrix(NA, nrow=length(groups_u), ncol=length(mod0$strata)))
    } else{
      max.times <- data.frame(X1=rep(NA, length(groups_u)))
    }

    # Find max times:
    for(group in 1:length(groups_u)){
      data_tmp <- data[data[,weight.names[1]]==groups_u[group],]

      mod_tmp <- rs.surv(formula, data_tmp, ratetable,
                         na.action, fin.date, method, conf.type,
                         conf.int, type, add.times, precision=1, rmap,
                         weight.table=-999) # add weights arg so that it is detected later for rmap

      mod_tmp <- summary(mod_tmp, times=sort(unique(mod0$time)))

      mod_l[[group]] <- mod_tmp

      if('strata' %in% names(mod0)){
        tajms_groupby <- aggregate(time~strata,
                                   data=data.frame(time=mod_tmp$time, strata=mod_tmp$strata),
                                   FUN = \(x) max(x))
        max.times[group,] <- tajms_groupby$time

      } else{
        max.times[group,] <- mod_tmp$time[length(mod_tmp$time)]
      }
    }

    # Find minimum among maximum follow-up times:
    find.max.time <- sapply(max.times, min)

    # Define data.frame:
    if('strata' %in% names(mod0)){
      find.max.time.df <- data.frame(final.time=find.max.time, strata=names(mod0$strata))
    } else{
      find.max.time.df <- data.frame(final.time=find.max.time)
    }

    # Find true length:
    if('strata' %in% names(mod0)){
      df.length <- data.frame(time=mod0$time, strata=rep(names(mod0$strata), mod0$strata))
      df.length <- merge(x=df.length, y=find.max.time.df, by='strata', all.x = TRUE)
    } else{
      df.length <- data.frame(time=mod0$time)
      df.length$final.time <- find.max.time.df$final.time
    }
    df.length <- df.length[df.length$time <= df.length$final.time,] # subset(df.length, time<=final.time)

    # Merge data sets:
    for(group in 1:length(groups_u)){
      mod_tmp <- mod_l[[group]]

      # Prepare data set for this group:
      if('strata' %in% names(mod0)){
        df_g <- data.frame(time=mod_tmp$time, surv=mod_tmp$surv, std.err=mod_tmp$std.err, strata=mod_tmp$strata)
        df_g <- merge(x=df_g, y=find.max.time.df, by='strata', all.x = TRUE)
      } else{
        df_g <- data.frame(time=mod_tmp$time, surv=mod_tmp$surv, std.err=mod_tmp$std.err)
        df_g$final.time <- find.max.time.df$final.time
      }

      df_g <- df_g[df_g$time <= df_g$final.time,] # subset(df_g, time<=final.time)

      # Take only relevant values:
      if('strata' %in% names(mod0)){
        df_g$st <- paste0(df_g$strata, df_g$time)
        df.length$st <- paste0(df.length$strata, df.length$time)
        df_g <- df_g[df_g$st %in% df.length$st,]
      } else{
        df_g <- df_g[df_g$time %in% df.length$time,]
      }

      if(group == 1){
        # Define objects:
        srvs <- matrix(NA, nrow=nrow(df.length), ncol=length(groups_u))
        std.errs <- matrix(NA, nrow=nrow(df.length), ncol=length(groups_u))
        # Save values:
        srvs[,group] <- df_g$surv*weight_val[group_val == groups_u[group]]
        std.errs[,group] <- df_g$std.err^2*weight_val[group_val == groups_u[group]]^2 # mod_tmp$surv[whi]^2*

      } else{
        # Save values:
        srvs[,group] <- df_g$surv*weight_val[group_val == groups_u[group]]
        std.errs[,group] <- df_g$std.err^2*weight_val[group_val == groups_u[group]]^2 # mod_tmp$surv[whi]^2*
      }
    }

    # Calculate the survival and std.err (this is sd(log(SNS))):
    out$surv <- rowSums(srvs)
    out$std.err <- sqrt(rowSums(std.errs)) / out$surv

    # Find which will stay in n.risk, n.event, n.censor:
    wh.stay <- which(out$time %in% df.length$time)

    # Take the correct times:
    out$time <- df.length$time

    # Calculate CI:
    se.fac <- sqrt(qchisq(conf.int, 1))						#factor needed for confidence interval
    if (conf.type == "plain") {
      out$lower <- as.vector(out$surv - out$std.err * se.fac * 		#surv + fac*se
                               out$surv)
      out$upper <- as.vector(out$surv + out$std.err * se.fac *
                               out$surv)
    }
    else if (conf.type == "log") {						#on log scale and back
      out$lower <- exp(as.vector(log(out$surv) - out$std.err *
                                   se.fac))
      out$upper <- exp(as.vector(log(out$surv) + out$std.err *
                                   se.fac))
    }
    else if (conf.type == "log-log") {						#on log-log scale and back
      out$lower <- exp(-exp(as.vector(log(-log(out$surv)) -
                                        out$std.err * se.fac/log(out$surv))))
      out$upper <- exp(-exp(as.vector(log(-log(out$surv)) +
                                        out$std.err * se.fac/log(out$surv))))
    }

    out$call <- call

    if('strata' %in% names(mod0)){
      out$strata <- table(df.length$strata)
    }
    # Remove unneeded objects:
    out$n.risk <- out$n.risk[wh.stay]
    out$n.event <- out$n.event[wh.stay]
    out$n.censor <- out$n.censor[wh.stay]

    return(out)

    # Regular part:
  } else{
    call <- match.call()
    if (!missing(rmap)) {
      if(is.null(weight.table)){
        rmap <- substitute(rmap)
      } else{
        if(weight.table != -999){
          rmap <- substitute(rmap)
        }
      }
    }
    rform <- rformulate(formula,data, ratetable, na.action,rmap)
    data <- rform$data								#the data set
    type <- match.arg(type, c("kaplan-meier", "fleming-harrington"))		#method of hazard -> survival scale transformation
    type <- match(type, c("kaplan-meier", "fleming-harrington"))
    method <- match.arg(method,c("pohar-perme", "ederer2", "hakulinen","ederer1"))	#method of relative surv. curve estimation
    method <- match(method,c("pohar-perme", "ederer2", "hakulinen","ederer1"))
    conf.type <- match.arg(conf.type,c("plain","log","log-log"))		#conf. interval type

    if (method == 3) {						#need potential follow-up time for Hak. method
      R <- rform$R
      coll <- match("year", attributes(ratetable)$dimid)
      year <- R[, coll]							#calendar year in the data
      if (missing(fin.date))
        fin.date <- max(rform$Y + year)					#final date for everybody set to the last day observed
      Y2 <- rform$Y							#change into potential follow-up time
      if (length(fin.date) == 1) 						#if final date equal for everyone
        Y2[rform$status == 1] <- fin.date - year[rform$status == 1]#set pot.time for those that died (equal to censoring time for others)
      else if (length(fin.date) == nrow(rform$R))
        Y2[rform$status == 1] <- fin.date[rform$status ==
                                            1] - year[rform$status == 1]
      else stop("fin.date must be either one value or a vector of the same length as the data")
      status2 <- rep(0, nrow(rform$X))						#stat2=0 for everyone
    }
    p <- rform$m								#number of covariates

    # strata - use same approach as in survival::survfit.formula:
    indx <- match(c("formula", "data"), names(call),
                  nomatch = 0)
    temp2 <- call[c(1, indx)]
    temp2[[1L]] <- quote(stats::model.frame)
    mf <- eval.parent(temp2)
    mf <- mf[rform$keep,]
    n <- nrow(mf)
    Terms2 <- terms(formula, c("strata", "cluster"))
    ll <- attr(Terms2, "term.labels")
    # Check if ratetable in linear predictor, if yes remove it:
    if(any(grepl(pattern = 'ratetable', x = ll))){
      ll <- ll[!grepl(pattern = 'ratetable', x = ll)]
    }
    if (length(ll) == 0)
      data$Xs <- factor(rep(1, n))
    else data$Xs <- strata(mf[ll])

    # OLD
    # if (p > 0) 									#if covariates
    #   data$Xs <- strata(rform$X[, ,drop=FALSE ])				#make strata according to covariates
    # else data$Xs <- rep(1, nrow(data))						#if no covariates, just put 1

    se.fac <- sqrt(qchisq(conf.int, 1))						#factor needed for confidence interval
    out <- NULL
    out$n <- table(data$Xs)							#table of strata
    out$time <- out$n.risk <- out$n.event <- out$n.censor <- out$surv <- out$std.err <- out$strata <-  NULL
    #out$index <- out$strata0 <- NULL
    # out$index = indices of the original times from the data among the times used for calculations
    # out$strata0 = the same as out$strata but only on the original times from the data


    for (kt in 1:length(out$n)) {						#for each stratum
      inx <- which(data$Xs == names(out$n)[kt])				#individuals within this stratum
      tis <- sort(unique(rform$Y[inx])) #unique times

      #if (method == 1 & all.times == TRUE) tis <- sort(union(rform$Y[inx],as.numeric(1:max(floor(rform$Y[inx])))))	#1-day long intervals used - to take into the account the continuity of the pop. part
      if (method == 1 & !missing(add.times)){
        #tis <- sort(union(rform$Y[inx],as.numeric(1:max(floor(rform$Y[inx])))))	#1-day long intervals used - to take into the account the continuity of the pop. part
        add.times <- pmin(as.numeric(add.times),max(rform$Y[inx]))
        tis <- sort(union(rform$Y[inx],as.numeric(add.times)))	#1-day long intervals used - to take into the account the continuity of the pop. part
      }
      if(method==3)tis <- sort(unique(pmin(max(tis),c(tis,Y2[inx]))))				#add potential times in case of Hakulinen
      #out$index <- c(out$index, which(tis %in% rform$Y[inx])+length(out$time))

      temp <- exp_prep(rform$R[inx,,drop=FALSE],rform$Y[inx],rform$ratetable,rform$status[inx],times=tis,fast=(method<3),prec=precision, ys=rform$start[inx])	#calculate the values for each interval of time

      out$time <- c(out$time, tis)						#add times
      out$n.risk <- c(out$n.risk, temp$yi)					#add number at risk for each time
      out$n.event <- c(out$n.event, temp$dni)					#add number of events for each time
      out$n.censor <- c(out$n.censor, temp$dci)
      # out$n.censor <- c(out$n.censor,  c(-diff(temp$yi),temp$yi[length(temp$yi)]) - temp$dni) 	#add number of censored for each time

      if(method==1){ 								#pohar perme method
        #approximate1 <- (temp$yidlisi/temp$yisi +temp$yidlisitt/temp$yisitt)/2
        #approximate <- (temp$yidlisiw/temp$yisi +temp$yidlisiw/temp$yisitt)/2		#approximation for integration
        approximate <- temp$yidlisiw
        #haz <- temp$dnisi/temp$yisi - temp$yidlisi/temp$yisi		#cumulative hazard increment on each interval
        haz <- temp$dnisi/temp$yisi - approximate			#cumulative hazard increment on each interval
        out$std.err <- c(out$std.err, sqrt(cumsum(temp$dnisisq/(temp$yisi)^2)))  #standard error on each interval
      }
      else if(method==2){							#ederer2 method
        haz <- temp$dni/temp$yi - temp$yidli/temp$yi			#cumulative hazard increment on each interval
        out$std.err <- c(out$std.err, sqrt(cumsum(temp$dni/(temp$yi)^2)))  #standard error on each interval
      }
      else if(method==3){							#Hakulinen method
        # Note: ys is not included in the following exp_prep; if Surv(start,stop,event) in formula it is not considered (you would have to use netweiDM probably):
        temp2 <- exp_prep(rform$R[inx,,drop=FALSE],Y2[inx],ratetable,status2[inx],times=tis)	#calculate the values for each interval of time
        popsur <- exp(-cumsum(temp2$yisidli/temp2$yisis))			#population survival
        haz <- temp$dni/temp$yi						#observed hazard on each interval
        out$std.err <- c(out$std.err, sqrt(cumsum(temp$dni/(temp$yi)^2)))  #standard error on each interval
      }
      else if(method==4){							#Ederer I
        popsur <- temp$sis/length(inx)					#population survival
        haz <- temp$dni/temp$yi						#observed hazard on each interval
        out$std.err <- c(out$std.err, sqrt(cumsum(temp$dni/(temp$yi)^2)))  #standard error on each interval
      }
      if(type==2)survtemp <- exp(-cumsum(haz))
      else survtemp <-  cumprod(1-haz)
      if(method>2){
        survtemp <- survtemp/popsur
      }
      out$surv <- c(out$surv,survtemp)
      out$strata <- c(out$strata, length(tis))				#number of times in this strata
      #out$strata0 <- c(out$strata0, length(unique(rform$Y[inx])))
    }
    if (conf.type == "plain") {
      out$lower <- as.vector(out$surv - out$std.err * se.fac * 		#surv + fac*se
                               out$surv)
      out$upper <- as.vector(out$surv + out$std.err * se.fac *
                               out$surv)
    }
    else if (conf.type == "log") {						#on log scale and back
      out$lower <- exp(as.vector(log(out$surv) - out$std.err *
                                   se.fac))
      out$upper <- exp(as.vector(log(out$surv) + out$std.err *
                                   se.fac))
    }
    else if (conf.type == "log-log") {						#on log-log scale and back
      out$lower <- exp(-exp(as.vector(log(-log(out$surv)) -
                                        out$std.err * se.fac/log(out$surv))))
      out$upper <- exp(-exp(as.vector(log(-log(out$surv)) +
                                        out$std.err * se.fac/log(out$surv))))
    }

    names(out$strata) <-  names(out$n)
    out$n <- as.vector(out$n)

    #names(out$strata0) <- names(out$n)

    if (p == 0){
      out$strata <-  NULL						#if no covariates
    }

    #if (method != 1) out$index <- out$strata0 <- NULL # if method != pohar-perme
    out$conf.type <- conf.type
    out$conf.int <- conf.int
    out$method <- method
    out$call <- call
    out$type <- rform$type
    class(out) <- c("survfit", "rs.surv")
    out
  }
}


#' Net Expected Sample Size Is Estimated
#'
#' Calculates how the sample size decreases in time due to population mortality
#'
#' The function calculates the sample size we can expect at a certain time
#' point if the patients die only due to population causes (population survival
#' * initial sample size in a certain category), i.e. the number of individuals
#' that remains at risk at given timepoints after the individuals who die due
#' to population causes are removed. The result should be used as a guideline
#' for the sensible length of follow-up interval when calculating the net
#' survival.
#'
#' The first column of the output reports the number of individuals at time 0.
#' The last column of the output reports the conditional expected (population)
#' survival time for each subgroup.
#'
#' @param formula a formula object, same as in \code{rs.surv}. The right-hand
#' side of the formula object includes the variable that defines the subgroups
#' (a variable of type \code{factor}) by which the expected sample size is to
#' be calculated.
#' @param data a data.frame in which to interpret the variables named in the
#' \code{formula}.
#' @param ratetable a table of event rates, organized as a \code{ratetable}
#' object, such as \code{slopop}.
#' @param times Times at which the calculation should be evaluated - in years!
#' @param rmap an optional list to be used if the variables are not organized
#' and named in the same way as in the \code{ratetable} object. See details of
#' the \code{rs.surv} function.
#' @return A list of values.
#' @seealso \code{rs.surv}
#' @references Pohar Perme, M., Pavlic, K. (2018) "Nonparametric Relative
#' Survival Analysis with the R Package relsurv". Journal of Statistical
#' Software. 87(8), 1-27, doi: "10.18637/jss.v087.i08"
#' @keywords survival
#' @examples
#'
#' data(slopop)
#' data(rdata)
#' rdata$agegr <-cut(rdata$age,seq(40,95,by=5))
#' nessie(Surv(time,cens)~agegr,rmap=list(age=age*365.241),
#' 	ratetable=slopop,data=rdata,times=c(1,3,5,10,15))
#'
nessie <- function (formula = formula(data), data = parent.frame(), ratetable = relsurv::slopop,times,rmap)

  #formula: for example Surv(time,cens)~sex
  #data: the observed data set
  #ratetable: the population mortality tables
  #times: the times at which to report NESS, if no default, then all unique times

{

  call <- match.call()
  if (!missing(rmap)) {
    rmap <- substitute(rmap)
  }
  na.action <- NA		#set the object just to be able to execute the rformulate call
  rform <- rformulate(formula, data, ratetable,na.action, rmap)			#get the data ready

  templab <- attr(rform$Terms,"term.labels")
  if(!is.null(attr(rform$Terms,"specials")$ratetable))templab <- templab[-length(templab)]	#delete the last term in the formula if the ratetable argument is there
  nameslist <- vector("list",length(templab))
  for(it in 1:length(nameslist)){
    valuetab <- table(data[,match(templab[it],names(data))])
    nameslist[[it]] <- paste(templab[it],names(valuetab),sep="")
  }
  names(nameslist) <- templab


  data <- rform$data								#the data set


  p <- rform$m								#number of covariates
  if (p > 0) 	{	#if covariates

    data$Xs <- my.strata(rform$X[,,drop=F],nameslist=nameslist)				#make strata according to covariates
    #data$Xs <- factor(data$Xs,levels=nameslist)				#order them in the same way as namelist
  }
  else data$Xs <- rep(1, nrow(data))						#if no covariates, just put 1

  if(!missing(times)) tis <- times
  else tis <- unique(sort(floor(rform$Y/365.241)))			#unique years of follow-up
  tis <- unique(c(0,tis))
  tisd <- tis*365.241


  out <- NULL
  out$n <- table(data$Xs)							#table of strata
  out$sp <- out$strata <- NULL
  # for (kt in 1:length(out$n)) {						#for each stratum
  for (kt in order(names(table(data$Xs)))) {						#for each stratum
    inx <- which(data$Xs == names(out$n)[kt])				#individuals within this stratum

    temp <- exp_prep(rform$R[inx,,drop=FALSE],rform$Y[inx],rform$ratetable,rform$status[inx],times=tisd,fast=FALSE)	#calculate the values for each interval of time

    out$time <- c(out$time, tisd)						#add times
    out$sp <- c(out$sp, temp$sis)						#add expected number of individuals alive
    out$strata <- c(out$strata, length(tis))				#number of times in this strata

    temp <- exp_prep(rform$R[inx,,drop=FALSE],rform$Y[inx],rform$ratetable,rform$status[inx],times=(seq(0,100,by=.5)*365.241)[-1],fast=FALSE)	#calculate the values for each interval of time
    out$povp <- c(out$povp,mean(temp$sit/365.241))
  }

  names(out$strata) <- names(out$n)[order(names(table(data$Xs)))]
  if (p == 0) out$strata <- NULL						#if no covariates

  mata <- matrix(out$sp,ncol=length(tis),byrow=TRUE)
  mata <- data.frame(mata)
  mata <- cbind(mata,out$povp)
  row.names(mata) <- names(out$n)[order(names(table(data$Xs)))]
  names(mata) <- c(tis,"c.exp.surv")


  cat("\n")
  print(round(mata,1))
  cat("\n")

  out$mata <- mata
  out$n <- as.vector(out$n)
  class(out) <- "nessie"
  invisible(out)
}







rs.period <- function (formula = formula(data), data = parent.frame(), ratetable = relsurv::slopop,
                       na.action, fin.date, method = "pohar-perme", conf.type = "log",
                       conf.int = 0.95,type="kaplan-meier",winst,winfin,diag.date,rmap)

  #formula: for example Surv(time,cens)~sex
  #data: the observed data set
  #ratetable: the population mortality tables
  #conf.type: confidence interval calculation (plain, log or log-log)
  #conf.int: confidence interval
  #winst: start of the period window (inclusive)
  #winfin: end of the period window (inclusive)

{

  call <- match.call()
  if (!missing(rmap)) {
    rmap <- substitute(rmap)
  }
  rform <- rformulate(formula, data, ratetable, na.action,rmap)			#get the data ready
  data <- rform$data								#the data set
  type <- match.arg(type, c("kaplan-meier", "fleming-harrington"))		#method of hazard -> survival scale transformation
  type <- match(type, c("kaplan-meier", "fleming-harrington"))
  method <- match.arg(method,c("pohar-perme", "ederer2", "hakulinen","ederer1"))	#method of relative surv. curve estimation
  method <- match(method,c("pohar-perme", "ederer2", "hakulinen","ederer1"))
  conf.type <- match.arg(conf.type,c("plain","log","log-log"))		#conf. interval type

  #machinations needed for period survival:
  R <- rform$R
  coll <- match("year", attributes(ratetable)$dimid)
  year <- R[, coll]							#calendar year in the data

  ys <- as.numeric(winst - year)
  yf <- as.numeric(winfin - year)

  relv <- which(ys <= rform$Y & yf>0)					#relevant individuals -> live up to the period window and were diagnosed before window end
  centhem <- which(yf < rform$Y)					#censor these - their event happens outside of the period window

  rform$status[centhem] <- 0
  rform$Y[centhem] <- yf[centhem]

  rform$Y <- rform$Y[relv]
  rform$X <- rform$X[relv,,drop=F]
  rform$R <- rform$R[relv,,drop=F]
  rform$status <- rform$status[relv]
  data <- data[relv,,drop=F]
  ys <- ys[relv]
  yf <- yf[relv]
  year <- year[relv]

  if (method == 3) {								#need potential follow-up time for Hak. method
    if (missing(fin.date))
      fin.date <- max(rform$Y + year)					#final date for everybody set to the last day observed
    Y2 <- rform$Y							#change into potential follow-up time
    if (length(fin.date) == 1) 						#if final date equal for everyone
      Y2[rform$status == 1] <- fin.date - year[rform$status == 1]#set pot.time for those that died (equal to censoring time for others)
    else if (length(fin.date[relv]) == nrow(rform$R)) 	{
      fin.date <- fin.date[relv]
      Y2[rform$status == 1] <- fin.date[rform$status ==
                                          1] - year[rform$status == 1]
    }
    else stop("fin.date must be either one value of a vector of the same length as the data")
    status2 <- rep(0, nrow(rform$X))						#stat2=0 for everyone
  }
  p <- rform$m								#number of covariates
  if (p > 0) 									#if covariates
    data$Xs <- strata(rform$X[, ,drop=FALSE ])				#make strata according to covariates
  else data$Xs <- rep(1, nrow(data))						#if no covariates, just put 1

  se.fac <- sqrt(qchisq(conf.int, 1))						#factor needed for confidence interval
  out <- NULL
  out$n <- table(data$Xs)							#table of strata
  out$time <- out$n.risk <- out$n.event <- out$n.censor <- out$surv <- out$std.err <- out$strata <- NULL
  for (kt in 1:length(out$n)) {						#for each stratum
    inx <- which(data$Xs == names(out$n)[kt])				#individuals within this stratum

    tis <- sort(unique(rform$Y[inx]))					#unique times
    if(method==3)tis <- sort(unique(pmin(max(tis),c(tis,Y2[inx]))))				#add potential times in case of Hakulinen

    ys <- pmax(ys,0)
    #tis <- sort(unique(c(tis,ys[ys>0]-1,ys[ys>0])))
    tis <- sort(unique(c(tis,ys[ys>0])))
    tis <- sort(unique(c(tis,tis-1,tis+1)))					#the day after exiting, the day before entering
    tis <- tis[-length(tis)]							#exclude the largest since it is beyond observation time (1 day later)

    temp <- exp_prep(rform$R[inx,,drop=FALSE],rform$Y[inx],rform$ratetable,rform$status[inx],times=tis,fast=(method<3),ys=ys)	#calculate the values for each interval of time

    out$time <- c(out$time, tis)						#add times
    out$n.risk <- c(out$n.risk, temp$yi)					#add number at risk for each time
    out$n.event <- c(out$n.event, temp$dni)					#add number of events for each time
    out$n.censor <- c(out$n.censor,  c(-diff(temp$yi),temp$yi[length(temp$yi)]) - temp$dni) 	#add number of censored for each time

    if(method==1){ 								#pohar perme method
      haz <- temp$dnisi/temp$yisi - temp$yidlisi/temp$yisi			#cumulative hazard increment on each interval
      out$std.err <- c(out$std.err, sqrt(cumsum(temp$dnisisq/(temp$yisi)^2)))  #standard error on each interval
    }
    else if(method==2){							#ederer2 method
      haz <- temp$dni/temp$yi - temp$yidli/temp$yi			#cumulative hazard increment on each interval
      out$std.err <- c(out$std.err, sqrt(cumsum(temp$dni/(temp$yi)^2)))  #standard error on each interval
    }
    else if(method==3){							#Hakulinen method
      temp2 <- exp_prep(rform$R[inx,,drop=FALSE],Y2[inx],rform$ratetable,status2[inx],times=tis,ys=ys)	#calculate the values for each interval of time
      popsur <- exp(-cumsum(temp2$yisidli/temp2$yisis))			#population survival
      haz <- temp$dni/temp$yi						#observed hazard on each interval
      out$std.err <- c(out$std.err, sqrt(cumsum(temp$dni/(temp$yi)^2)))  #standard error on each interval
    }
    else if(method==4){							#Ederer I
      popsur <- temp$sis/length(inx)					#population survival
      haz <- temp$dni/temp$yi						#observed hazard on each interval
      out$std.err <- c(out$std.err, sqrt(cumsum(temp$dni/(temp$yi)^2)))  #standard error on each interval
    }
    if(type==2)survtemp <- exp(-cumsum(haz))
    else survtemp <-  cumprod(1-haz)
    if(method>2){
      survtemp <- survtemp/popsur
    }
    out$surv <- c(out$surv,survtemp)
    out$strata <- c(out$strata, length(tis))				#number of times in this strata
  }
  if (conf.type == "plain") {
    out$lower <- as.vector(out$surv - out$std.err * se.fac * 		#surv + fac*se
                             out$surv)
    out$upper <- as.vector(out$surv + out$std.err * se.fac *
                             out$surv)
  }
  else if (conf.type == "log") {						#on log scale and back
    out$lower <- exp(as.vector(log(out$surv) - out$std.err *
                                 se.fac))
    out$upper <- exp(as.vector(log(out$surv) + out$std.err *
                                 se.fac))
  }
  else if (conf.type == "log-log") {						#on log-log scale and back
    out$lower <- exp(-exp(as.vector(log(-log(out$surv)) -
                                      out$std.err * se.fac/log(out$surv))))
    out$upper <- exp(-exp(as.vector(log(-log(out$surv)) +
                                      out$std.err * se.fac/log(out$surv))))
  }
  names(out$strata) <- names(out$n)
  if (p == 0) out$strata <- NULL						#if no covariates
  out$n <- as.vector(out$n)
  out$conf.type <- conf.type
  out$conf.int <- conf.int
  out$method <- method
  out$call <- call
  out$type <- rform$type
  class(out) <- c("survfit", "rs.surv")
  out
}



#' expprep2 function
#'
#' Helper calculation function using C code. Saved also as exp_prep (unexported
#' function).
#'
#' Helper function used in rs.surv and other relsurv functions.
#'
#' @param x matrix of demographic covariates - each individual has one line
#' @param y follow-up time for each individual (same length as nrow(x))
#' @param ratetable rate table used for calculation
#' @param status status for each individual (same length as nrow(x)!), not
#' needed if we only need Spi, status needed for rs.surv
#' @param times times at which we wish to evaluate the quantities, not needed
#' if we only need Spi, times needed for rs.surv
#' @param fast for mpp method only
#' @param ys entry times (if empty, individuals are followed from time 0)
#' @param prec deprecated
#' @param cmp should cmpfast.C be used
#' @param netweiDM should new netwei script be used
#' @return List containing the calculated hazards and probabilities using the
#' population mortality tables.
#' @seealso rs.surv
#' @keywords survival
#' @export expprep2
expprep2 <- function (x, y,ratetable,status,times,fast=FALSE,ys,prec,cmp=F,netweiDM=FALSE) {			#function that prepares the data for C function call

  #x= matrix of demographic covariates - each individual has one line
  #y= follow-up time for each individual (same length as nrow(x)!)
  #ratetable= rate table used for calculation
  #status= status for each individual (same length as nrow(x)!), not needed if we only need Spi, status needed for rs.surv
  #times= times at which we wish to evaluate the quantities, not needed if we only need Spi, times needed for rs.surv
  #fast=for mpp method only
  #netweiDM=should new netwei script be used

  x <- as.matrix(x)
  if (ncol(x) != length(dim(ratetable)))
    stop("x matrix does not match the rate table")
  atts <- attributes(ratetable)

  cuts <- atts$cutpoints

  if (is.null(atts$type)) {
    rfac <- atts$factor
    us.special <- (rfac > 1)
  }
  else {
    rfac <- 1 * (atts$type == 1)
    us.special <- (atts$type == 4)
  }
  if (length(rfac) != ncol(x))
    stop("Wrong length for rfac")


  if (any(us.special)) {
    if (sum(us.special) > 1)
      stop("Two columns marked for special handling as a US rate table")
    cols <- match(c("age", "year"), atts$dimid)
    if (any(is.na(cols)))
      stop("Ratetable does not have expected shape")
    if (exists("as.Date")) {
      bdate <- as.Date("1960/1/1") + (x[, cols[2]] - x[,
                                                       cols[1]])
      byear <- format(bdate, "%Y")
      offset <- as.numeric(bdate - as.Date(paste(byear,
                                                 "01/01", sep = "/")))
    } else stop("Can't find an appropriate date class\n")
    # Tole je bilo prej:
    # else if (exists("date.mdy")) {
    #   bdate <- as.date(x[, cols[2]] - x[, cols[1]])
    #   byear <- date.mdy(bdate)$year
    #   offset <- bdate - mdy.date(1, 1, byear)
    # }
    x[, cols[2]] <- x[, cols[2]] - offset
    if (any(rfac > 1)) {
      temp <- which(us.special)
      nyear <- length(cuts[[temp]])
      nint <- rfac[temp]
      cuts[[temp]] <- round(approx(nint * (1:nyear), cuts[[temp]],
                                   nint:(nint * nyear))$y - 1e-04)
    }
  }

  if(!missing(status)){		#the function was called from rs.surv
    if(length(status)!=nrow(x))    stop("Wrong length for status")

    if(missing(times))    times <- sort(unique(y))

    if (any(times < 0))
      stop("Negative time point requested")
    ntime <- length(times)
    if(missing(ys)) ys <- rep(0,length(y))
    #    times2 <- times
    #    times2[1] <- preci

    # It may be an integer...check and make it numeric. Otherwise problems in C functions:
    if(is.integer(x[1,1])){
      x <- apply(x, 2, as.numeric)
    }

    if(cmp)   temp <- .Call("cmpfast",  as.integer(rfac), 		#fast=pohar-perme or ederer2 - data from pop. tables only while under follow-up
                            as.integer(atts$dim), as.double(unlist(cuts)), ratetable,
                            x, y, ys,as.integer(status), times,PACKAGE="relsurv")
    else if(fast&!missing(prec))    temp <- .Call("netfastpinter2",  as.integer(rfac), 		#fast=pohar-perme or ederer2 - data from pop. tables only while under follow-up
                                                  as.integer(atts$dim), as.double(unlist(cuts)), ratetable,
                                                  x, y, ys,as.integer(status), times,prec,PACKAGE="relsurv")
    else if(fast&missing(prec))    temp <- .Call("netfastpinter",  as.integer(rfac), 		#fast=pohar-perme or ederer2 - data from pop. tables only while under follow-up
                                                 as.integer(atts$dim), as.double(unlist(cuts)), ratetable,
                                                 x, y, ys,as.integer(status), times,PACKAGE="relsurv")
    else if(netweiDM==TRUE)    temp <- .Call("netweiDM",  as.integer(rfac),
                                             as.integer(atts$dim), as.double(unlist(cuts)), ratetable,
                                             x, y, ys,as.integer(status), times,PACKAGE="relsurv")
    else    temp <- .Call("netwei",  as.integer(rfac),
                          as.integer(atts$dim), as.double(unlist(cuts)), ratetable,
                          x, y, as.integer(status), times,PACKAGE="relsurv")
  }
  else{				#only expected survival at time y is needed for each individual
    if(length(y)==1)y <- rep(y,nrow(x))
    if(length(y)!=nrow(x)) stop("Wrong length for status")
    temp <- .Call("expc",  as.integer(rfac),
                  as.integer(atts$dim), as.double(unlist(cuts)), ratetable,
                  x, y,PACKAGE="relsurv")
    temp  <- temp$surv
  }
  temp
}
