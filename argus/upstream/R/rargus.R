
logpsi <- function(chi){ pgamma(0.5*chi^2,shape=1.5,log.p=TRUE)+log(0.5)}
dargus <- function(x, chi, log = FALSE) {
  ## ------------------------------------------------------------------------
  ## return density of argus distributed random variates
  ##
  ## density proportional to
  ## 
  ##      1>= x >= 0
  ## ------------------------------------------------------------------------
  ## Arguments:
  ##
  ##   x ....... values to be evaluated
  ##   chi   ... parameter for distribution                                   
  ##   log   ... if TRUE the logarithm of the density will be returned
  ## ------------------------------------------------------------------------
        y <- 1.0 - x*x
        A = 3*log(chi) - log(sqrt(2*pi)) - logpsi(chi)
        if(log) return( A + log(x) + 0.5*log1p(-x*x) - chi^2 * (1-x*x) * 0.5 )
        return( exp( A + log(x) + 0.5*log1p(-x*x) - chi^2 * (1-x*x) * 0.5 ) )
}




pargus <-
function(x, chi, lower=TRUE, log.p = FALSE) {
  ## ------------------------------------------------------------------------
  ## return CDF of argus distributed random variates
  ##
  ## density proportional to
  ##    TODO add density formula
  ##      1>= x >= 0
  ## ------------------------------------------------------------------------
  ## Arguments:
  ##
  ##   x ....... values tobe evaluated
  ##   chi   ... parameter for distribution                                   
  ##   lower ... if FALSE 1-CDF is returned
  ##   log   ... if TRUE the logarithm of the CDF will be returned
  ## ------------------------------------------------------------------------
  ## below code based on formula 
  ## CDF(x,chi) = 1.-pgamma(0.5*chi*chi*(1-x*x),shape=1.5)/pgamma(0.5*chi*chi,shape=1.5)
  reslogupper <-   pgamma(0.5*chi*chi*(1-x*x),shape=1.5,log.p=TRUE)-pgamma(0.5*chi*chi,shape=1.5,log.p=TRUE)
  if(!lower & log.p) return(reslogupper) 
  if(!lower & !log.p) p <- return(exp(reslogupper))
  if(!log.p) return(-expm1(reslogupper)) # case lower and !log.p
  return(log(-expm1(reslogupper))) # case lower and log.p
}

#######################
qargus <-
function(p, chi, lower=TRUE, log.p = FALSE) {
  ## ------------------------------------------------------------------------
  ## return CDF of argus distributed random variates
  ##
  ## density proportional to
  ##    TODO add density formula
  ##      1>= x >= 0
  ## ------------------------------------------------------------------------
  ## Arguments:
  ##
  ##   p ....... probabilities for which the quantiles are to be evaluated
  ##   chi   ... parameter for distribution                                   
  ##   lower ... if FALSE F^(-1)(1-p) is returned
  ##   log   ... if TRUE the logarithm of the probability will be returned
  ## ------------------------------------------------------------------------

  ## generate sample
  
  if(log.p) p <- exp(p)
  C <- pgamma(chi*chi*0.5,shape=1.5)
  if(lower) return(sqrt(1-qgamma((1-p)*C,shape=1.5)*2/chi^2) )
  return(sqrt(1-qgamma(p*C,shape=1.5)*2/chi^2) )
}


###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################

#rargus <- function (n = length(chi), chi, method = c("RoU", "inversion")) 
#{
#    chi <- rep(chi, n%/%length(chi))
#    if (method[1] == "RoU") 
#        return(.Call("rargus", n = n%/%length(chi), chi))
#    return(rargus.inversion(chi))
#}



rargus <-
function(n=length(chi), chi, method=c("inversion","RoU")) {
  ## ------------------------------------------------------------------------
  ## Generate ARGUS distributed random variates
  ##
  ## density proportional to
  ## 
  ##       1 >= x >= 0
  ## ------------------------------------------------------------------------
  ## Arguments:
  ##
  ##   n ....... total sample size, must be a multiple of length chi
  ##   chi   ... parameter for distribution                                   
  ## ------------------------------------------------------------------------

  ## generate sample
  if(method[1]=="inversion"){
    if(length(chi)==1 ) return( rargus.inversion.1chi(n,chi) )
    chi <- rep(chi,n%/%length(chi))
    return( rargus.inversion(chi) )
  }
  chi <- rep(chi,n%/%length(chi))
  return( .Call("rargusRoU", n=n%/%length(chi), chi) )# calls RoU method
}
#rargusN(n=1, chi=0.5, method=c("fastest")) 




rargus.inversion.1chi <- function(n,chi){
# generates n argus Random variate for the single value chi
# chi ... vector of the parameters chi of the argus distribution

if( chi>=1){
#  if(!exists("argus.env$_argus.genora")) assign("_argus.genora", rargus.setup(uresolution=1.e-13), envir = argus.env) #.GlobalEnv
#  genora <- get("_argus.genora",envir=argus.env)
  CU <- argus.env$genora$splG1to10(chi)*(1-runif(n)) 
  return( sqrt(1-2* Runuran::uq(argus.env$genora$gH0,CU) /(chi*chi)))
}else if( chi>=.1){
#  if(!exists("argus.env$_argus.genora")) assign("_argus.genora", rargus.setup(uresolution=1.e-13), envir = argus.env) 
#  genora <- get("_argus.genora",envir=argus.env)
  CU <- argus.env$genora$splG01(chi)*(1-runif(n)) 
  return( sqrt(1-2* Runuran::uq(argus.env$genora$gH1,CU/argus.env$genora$C1) /(chi*chi)))
}else if( chi>=.01){
#  if(!exists("argus.env$_argus.genora")) assign("_argus.genora", rargus.setup(uresolution=1.e-13), envir = argus.env)
#  genora <- get("_argus.genora",envir=argus.env)
  CU <- argus.env$genora$splG01(chi)*(1-runif(n)) 
  return( sqrt(1-2* Runuran::uq(argus.env$genora$gH2,CU/argus.env$genora$C2) /(chi*chi)))
}
Y <- 0.5*chi^2*(1-runif(n))^(2/3)
if(chi > 1.e-5){ 
 Y<- ( (Y*(1+Y*(1+Y*(0.5+Y/6))))*(1/3-0.1*Y+2*(1/3-0.1*chi^2))  -0.5*Y*Y*(1+Y/3) ) } # uses approximation (16) of the paper

  return( sqrt(1-2*Y /(chi*chi)) )

#
#
#  gen <- Runuran::pinv.new(pdf=function(x)log(x)+0.5*log(1-x*x)-0.5*chi*chi*(1-x*x),lb=0,ub=1,islog=T,uresolution=1.e-12)
#  return(Runuran::uq(gen,runif(n)))
}
#system.time( for(i in 1:1000){ y<-rargus.inversion.1chi(1.e3,chi=1) } )



rargus.inversion <- function(chi){
# generates one argus Random variate for each of the values in chi
# chi ... vector of the parameters chi of the argus distribution

#if(!exists("argus.env$_argus.genora")) assign("_argus.genora", rargus.setup(uresolution=1.e-13), envir = argus.env)
#genora <- get("_argus.genora",envir=argus.env)

rgam1.5chi.gr.1 <- function(chi){
#for chi values >1
CU <- argus.env$genora$splG1to10(chi)*(1-runif(length(chi))) 
return(Runuran::uq(argus.env$genora$gH0,CU))
}
rgam1.5chi..1to1 <- function(chi){
#for chi values in(0.1,1)
CU <- argus.env$genora$splG01(chi)*(1-runif(length(chi))) 
return(Runuran::uq(argus.env$genora$gH1,CU/argus.env$genora$C1))
}
rgam1.5chi..01to.1 <- function(chi){
#for chi values in(0.01,.1)
CU <- argus.env$genora$splG01(chi)*(1-runif(length(chi))) 
return(Runuran::uq(argus.env$genora$gH2,CU/argus.env$genora$C2))
}
rgam1.5chi.less.01 <- function(chi){
#for chi values lessthan 0.01
  YY <- 0.5*chi^2*(1-runif(length(chi)))^(2/3)
  iv <- chi>1.e-5
  Y <- YY[iv]
  YY[iv] <- ( (Y*(1+Y*(1+Y*(0.5+Y/6))))*(1/3-0.1*Y+2*(1/3-0.1*chi[iv]^2))  -0.5*Y*Y*(1+Y/3) ) # uses approximation (16) of the paper
  return(YY) }
#   return(ifelse(chi < 1.e-5,Y ,
#          (Y*(1+Y*(1+Y*(0.5+Y/6))))*(1/3-0.1*Y+sqrt(2*pi)*(1/3-0.1*chi^2)*sqrt(2/pi))  -0.5*Y*Y*(1+Y/3) ) ) # uses approximation (16) of the paper
##          (Y*(1+Y*(1+Y*(0.5+Y/6))))*(1/3-0.1*Y+sqrt(2*pi)*pgamma(0.5*(chi*chi),1.5)/chi^3)-0.5*Y*Y*(1+Y/3) ) )
Y <- numeric(length(chi))
lv <- chi>=1 #which(chi>=1)
Y[lv]<- rgam1.5chi.gr.1(chi[lv])
lv <- !chi>=1 & chi>=0.1 
Y[lv]<-  rgam1.5chi..1to1(chi[lv])
lv <- chi<0.1 &chi>=0.01
Y[lv]<-  rgam1.5chi..01to.1(chi[lv])
lv <- chi<0.01
Y[lv]<-  rgam1.5chi.less.01(chi[lv])
return(sqrt(1-2*Y/(chi*chi)))
}

rargus.setup <- function(uresolution=1.e-10){
# setup that returns the object required for generating argus by inversion
# has to be be run only once and works for any chi value 
# uresolution ... must not be smaller than 1.e-12
G <- function(x){pgamma(x,1.5)}
xv<- seq(0,1.1,0.001)
splG01 <-splinefun(xv,pgamma(0.5*(xv*xv),1.5))
xv<- seq(0.8,10,0.01)
splG1to10 <-splinefun(xv,pgamma(0.5*(xv*xv),1.5))
res <- list(
  gH0 = Runuran::pinv.new(pdf=function(x)0.5*log(x)-x,lb=0,ub=Inf,islog=T,center=0.5,ures=uresolution),
  gH1 = Runuran::pinv.new(pdf=function(x)0.5*log(x)-x,lb=0,ub=0.5,islog=T,ures=uresolution*0.001),
  gH2 = Runuran::pinv.new(pdf=function(x)0.5*log(x)-x,lb=0,ub=0.005,center=0.004,islog=T,ures=uresolution*0.001),
  C1 = G(0.5),
  C2 = G(0.005),
  splG01 = splG01,
  splG1to10 = splG1to10  )
return(res)
}

# environment to store the tables necessary for fast inversion
argus.env <- new.env(parent = emptyenv())

.onAttach <- function(libname,pkggname) {  
# calcuates the tables necessary for rargus.inversion()
# and stores them in the object argus.env$genora
  assign("genora", rargus.setup(uresolution=1.e-10), envir = argus.env)
}




