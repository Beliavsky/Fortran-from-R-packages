#'Finds the intersection of two trawl sets
#'@name fit_trawl_intersection
#'@param fct1 specifies the type of the first trawl function
#'@param fct2 specifies the type of the second trawl function
#'@param lambda11,lambda12,w1 parameters of the (double) exponential trawl
#'  functions of the first process
#'@param delta1,gamma1 parameters of the supIG trawl functions of the first
#'  process
#'@param alpha1,H1 parameters of the long memory trawl function of the first
#'  process
#'@param lambda21,lambda22,w2 parameters of the (double) exponential trawl
#'  functions of the second process
#'@param delta2,gamma2 parameters of the supIG trawl functions of the second
#'  process
#'@param alpha2,H2 parameters of the long memory trawl function of the second
#'  process
#'@param LM1 Lebesgue measure of the first trawl
#'@param LM2 Lebesgue measure of the second trawl
#'@param plotdiag binary variable specifying whether or not diagnostic plots
#'  should be provided
#'@return The Lebesgue measure of the intersection of the two trawl sets
#'@details Computes \eqn{R_{12}(0)=\mbox{Leb}(A_1 \cap A_2)} based on two trawl
#'  functions \eqn{g_1} and \eqn{g_2}.
#'@export

fit_trawl_intersection <- function(fct1 =base::c("Exp", "DExp", "supIG", "LM"),fct2 =base::c("Exp", "DExp", "supIG", "LM"),
                                   lambda11=0,lambda12=0,w1=0,delta1=0,gamma1=0,alpha1=0,H1=0,
                                   lambda21=0,lambda22=0,w2=0,delta2=0,gamma2=0,alpha2=0,H2=0,LM1, LM2,plotdiag=FALSE){


  if(fct1=="Exp"){    g1 <- function(z) {trawl_Exp(z,lambda11)}}
  if(fct1=="DExp"){   g1 <- function(z){trawl_DExp(z,w1,lambda11,lambda12)}}
  if(fct1 =="supIG"){ g1 <- function(z){trawl_supIG(z,delta1,gamma1)}}
  if(fct1 == "LM"){   g1 <- function(z){trawl_LM(z,alpha1, H1)}}

  if(fct2=="Exp"){    g2 <- function(z) {trawl_Exp(z,lambda21)}}
  if(fct2=="DExp"){   g2 <- function(z){trawl_DExp(z,w2,lambda21,lambda22)}}
  if(fct2 =="supIG"){ g2 <- function(z){trawl_supIG(z,delta2,gamma2)}}
  if(fct2 == "LM"){   g2 <- function(z){trawl_LM(z,alpha2,H2)}}

  g  <- function(z) {g1(z)-g2(z)}
  #Find the root(s)
  if(plotdiag){
    graphics::curve(g,-1000,0)
    graphics::plot((-1000:0),g((-1000:0)))
  }


  #Find all roots
  lb <- -1000
  up <- -0.000001
  allroots <- rootSolve::uniroot.all(g, base::c(lb, up))
  l <- base::length(allroots)
  if(l==0){#no roots
    base::print("There are no roots")
    #omax <- stats::optimize(g, interval=base::c(-500, 0), maximum=TRUE)
    #omin <- stats::optimize(g, interval=base::c(-500, 0), maximum=FALSE)

    #if((omax$objective>=0)&&(omin$objective>=0)){R12 <-LM2}
    #if((omax$objective<=0)&&(omin$objective<=0)){R12 <-LM1}
    R12<-base::min(LM1,LM2)
  }
  if(l==1){ #one root
    base::print("There is one root")

    if(g((lb+allroots)/2)<=0){
      R12 <- stats::integrate(g1,-Inf,allroots,stop.on.error=FALSE)$value +stats::integrate(g2,allroots,0,stop.on.error=FALSE)$value
    }
    if(g((lb+allroots)/2)>0){
      R12 <- stats::integrate(g2,-Inf,allroots,stop.on.error=FALSE)$value +stats::integrate(g1,allroots,0,stop.on.error=FALSE)$value
    }
  }
  if(l==2){ #two roots
    base::print("There are two roots")

    if(g((lb+allroots[1])/2)<=0){
      R12 <- stats::integrate(g1,-Inf,allroots[1],stop.on.error=FALSE)$value +stats::integrate(g2,allroots[1],allroots[2],stop.on.error=FALSE)$value + stats::integrate(g1,allroots[2],0,stop.on.error=FALSE)$value
    }
    if(g((lb+allroots[1])/2)>0){
      R12 <- stats::integrate(g2,-Inf,allroots[1],stop.on.error=FALSE)$value +stats::integrate(g1,allroots[1],allroots[2],stop.on.error=FALSE)$value + stats::integrate(g2,allroots[2],0,stop.on.error=FALSE)$value
    }
  }
  if(l>=3){#more roots?!?!!?!
    base::print("There are more than three roots")
    R12 <- base::min(LM1,LM2)

  }

  return(R12)
}


#'Finds the intersection of two long memory (LM) trawl sets
#'@name fit_trawl_intersection_LM
#'@param alpha1,H1,alpha2,H2 parameters of the two long memory trawls
#'@param LM1 Lebesgue measure of the first trawl
#'@param LM2 Lebesgue measure of the second trawl
#'@param plotdiag binary variable specifying whether or not diagnostic plots
#'  should be provided
#'@return the Lebesgue measure of the intersection of the two trawl sets
#'@details Computes \eqn{R_{12}(0)=\mbox{Leb}(A_1 \cap A_2)} based on two trawl
#'  functions \eqn{g_1} and \eqn{g_2}.
#'@export

fit_trawl_intersection_LM <- function(alpha1,H1,alpha2,H2,LM1, LM2,plotdiag=FALSE){

g1 <- function(z){trawl_LM(z,alpha1,H1)}
g2 <- function(z){trawl_LM(z,alpha2,H2)}

g  <- function(z) {g1(z)-g2(z)}
#Find the root(s)
if(plotdiag){
  graphics::curve(g,-1000,0)
  graphics::plot((-1000:0),g((-1000:0)))
}


#Find all roots
lb <- -1000
up <- -0.000001
allroots <- rootSolve::uniroot.all(g, base::c(lb, up))
l <- base::length(allroots)
if(l==0){#no roots
  base::print("There are no roots")
  omax <- stats::optimize(g, interval=base::c(-500, 0), maximum=TRUE)
  omin <- stats::optimize(g, interval=base::c(-500, 0), maximum=FALSE)

  if((omax$objective>=0)&&(omin$objective>=0)){R12 <-LM2}
  if((omax$objective<=0)&&(omin$objective<=0)){R12 <-LM1}
}
if(l==1){ #one root
  base::print("There is one root")

  if(g((lb+allroots)/2)<=0){
    R12 <- alpha1/(H1-1)*(1+(-allroots)/alpha1)^(1-H1) + alpha2/(H2-1)*(1-(1+(-allroots)/alpha2)^(1-H2))
    #integrate(g1,-Inf,allroots)$value +integrate(g2,allroots,0)$value
  }
  if(g((lb+allroots)/2)>0){
    #R12 <- integrate(g2,-Inf,allroots)$value +integrate(g1,allroots,0)$value
    R12 <- alpha2/(H2-1)*(1+(-allroots)/alpha2)^(1-H2) + alpha1/(H1-1)*(1-(1+(-allroots)/alpha1)^(1-H1))
  }
}
if(l==2){ #two roots
  base::print("There are two roots")

  if(g((lb+allroots[1])/2)<=0){
    R12 <- stats::integrate(g1,-Inf,allroots[1])$value +stats::integrate(g2,allroots[1],allroots[2])$value + stats::integrate(g1,allroots[2],0)$value
  }
  if(g((lb+allroots[1])/2)>0){
    R12 <- stats::integrate(g2,-Inf,allroots[1])$value +stats::integrate(g1,allroots[1],allroots[2])$value + stats::integrate(g2,allroots[2],0)$value
  }
}
if(l>=3){#more roots?!?!!?!
  print("There are more than three roots")
  R12 <- base::min(LM1,LM2)

}

return(R12)
}

#'Finds the intersection of two exponential trawl sets
#'@name fit_trawl_intersection_Exp
#'@param lambda1,lambda2 parameters of the two exponential trawls
#'@param LM1 Lebesgue measure of the first trawl
#'@param LM2 Lebesgue measure of the second trawl
#'@param plotdiag binary variable specifying whether or not diagnostic plots
#'  should be provided
#'@return The Lebesgue measure of the intersection of the two trawl sets
#'@details Computes \eqn{R_{12}(0)=\mbox{Leb}(A_1 \cap A_2)} based on two trawl
#'  functions \eqn{g_1} and \eqn{g_2}.
#'@export

fit_trawl_intersection_Exp <- function(lambda1,lambda2,LM1, LM2,plotdiag=FALSE){

  g1 <- function(z) {trawl_Exp(z,lambda1)}
  g2 <- function(z) {trawl_Exp(z,lambda2)}

  g  <- function(z) {g1(z)-g2(z)}
  #Find the root(s)
  if(plotdiag){
    graphics::curve(g,-1000,0)
    graphics::plot((-1000:0),g((-1000:0)))
  }


  #Find all roots
  lb <- -1000
  up <- -0.000001
  allroots <- rootSolve::uniroot.all(g, base::c(lb, up))
  l <- length(allroots)
  if(l==0){#no roots
    base::print("There are no roots")
    omax <- stats::optimize(g, interval=base::c(-500, 0), maximum=TRUE)
    omin <- stats::optimize(g, interval=base::c(-500, 0), maximum=FALSE)

    if((omax$objective>=0)&&(omin$objective>=0)){R12 <-LM2}
    if((omax$objective<=0)&&(omin$objective<=0)){R12 <-LM1}
  }
  if(l==1){ #one root
    base::print("There is one root")

    if(g((lb+allroots)/2)<=0){
      R12 <- stats::integrate(g1,-Inf,allroots)$value +stats::integrate(g2,allroots,0)$value
    }
    if(g((lb+allroots)/2)>0){
      R12 <- stats::integrate(g2,-Inf,allroots)$value +stats::integrate(g1,allroots,0)$value
    }
  }
  if(l==2){ #two roots
    print("There are two roots")

    if(g((lb+allroots[1])/2)<=0){
      R12 <- stats::integrate(g1,-Inf,allroots[1])$value +stats::integrate(g2,allroots[1],allroots[2])$value + stats::integrate(g1,allroots[2],0)$value
    }
    if(g((lb+allroots[1])/2)>0){
      R12 <- stats::integrate(g2,-Inf,allroots[1])$value +stats::integrate(g1,allroots[1],allroots[2])$value + stats::integrate(g2,allroots[2],0)$value
    }
  }
  if(l>=3){#more roots?!?!!?!
    base::print("There are more than three roots")
    R12 <- base::min(LM1,LM2)

  }

  return(R12)
}


