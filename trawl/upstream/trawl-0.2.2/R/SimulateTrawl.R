#'Evaluates the exponential trawl function
#'@name trawl_Exp
#'@param x the argument at which the exponential trawl function will be
#'  evaluated
#'@param lambda the parameter \eqn{\lambda} in the exponential trawl
#'@return The exponential trawl function evaluated at x
#'@details The trawl function is parametrised by parameter \eqn{\lambda > 0}  as
#'  follows: \deqn{g(x) = e^{\lambda x},  \mbox{ for }  x \le 0.}
#'@export
trawl_Exp <- function(x,lambda){exp(x*lambda)}

#'Evaluates the double exponential trawl function
#'@name trawl_DExp
#'@param x the argument at which the double exponential trawl function will be
#'  evaluated
#'@param w parameter in the double exponential trawl
#'@param lambda1 the parameter \eqn{\lambda_1} in the double exponential trawl
#'@param lambda2 the parameter \eqn{\lambda_2} in the double exponential trawl
#'@return The double exponential trawl function evaluated at x
#'@details The trawl function is parametrised by parameters \eqn{0\leq w\leq 1}
#'  and \eqn{\lambda_1, \lambda_2 > 0}  as follows: \deqn{g(x) = w e^{\lambda_1
#'  x}+(1-w) e^{\lambda_2 xz},  \mbox{ for }  x \le 0.}
#'@export
trawl_DExp <- function(x,w,lambda1,lambda2){w*exp(lambda1*x)+(1-w)*exp(lambda2*x)}


#'Evaluates the supIG trawl function
#'@name trawl_supIG
#'@param x the argument at which the supIG trawl function will be evaluated
#'@param delta the parameter \eqn{\delta} in the supIG trawl
#'@param gamma the parameter \eqn{\gamma} in the supIG trawl
#'@return The supIG trawl function evaluated at x
#'@details The trawl function is parametrised by the two parameters \eqn{\delta
#'  \geq 0} and \eqn{\gamma \geq 0} as follows: \deqn{gd(x) =
#'  (1-2x\gamma^{-2})^{-1/2}\exp(\delta \gamma(1-(1-2x\gamma^{-2})^{1/2})),
#'  \mbox{ for } x \le 0.} It is assumed that \eqn{\delta} and \eqn{\gamma} are
#'  not simultaneously equal to zero.
#'@export
trawl_supIG <-function(x,delta,gamma){
  r<-(1-2*x/gamma^2)^(-1/2)*exp(delta*gamma*(1-(1-2*x/gamma^2)^(1/2)))
  return(r)
}


#'Evaluates the long memory trawl function
#'@name trawl_LM
#'@param x the argument at which the long memory trawl function will be
#'  evaluated
#'@param alpha the parameter \eqn{\alpha} in the long memory trawl
#'@param H the parameter \eqn{H} in the long memory trawl
#'@return the long memory trawl function evaluated at x
#'@details The trawl function is parametrised by the two parameters \eqn{H> 1}
#'  and \eqn{\alpha > 0} as follows: \deqn{g(x) = (1-x/\alpha)^{-H},  \mbox{ for
#'  }  x \le 0.}
#'@export
trawl_LM <- function(x,alpha, H){(1-x/alpha)^(-H)}


#'Simulates a univariate trawl process
#'@name sim_UnivariateTrawl
#'@param t parameter which specifying the length of the time interval
#'  \eqn{[0,t]} for which a simulation of the trawl process is required.
#'@param Delta parameter \eqn{\Delta} specifying the length of the time step,
#'  the default is 1
#'@param burnin parameter specifying the length of the burn-in period at the
#'  beginning of the simulation
#'@param marginal parameter specifying the marginal distribution of the trawl
#'@param trawl parameter specifying the type of trawl function
#'@param v parameter of the Poisson distribution
#'@param m parameter of the negative binomial distribution
#'@param theta parameter \eqn{\theta} of the negative binomial distribution
#'@param lambda1 parameter \eqn{\lambda_1} of the exponential (or
#'  double-exponential) trawl function
#'@param lambda2 parameter \eqn{\lambda_2} of the double-exponential trawl
#'  function
#'@param w parameter of the double-exponential trawl function
#'@param delta parameter \eqn{\delta} of the supIG trawl function
#'@param gamma parameter \eqn{\gamma} of the supIG trawl function
#'@param alpha parameter \eqn{\alpha} of the long memory trawl function
#'@param H parameter of the long memory trawl function
#'@details This function simulates a univariate trawl process with either
#'  Poisson or negative binomial marginal law. For the trawl function there are
#'  currently four choices: exponential, double-exponential, supIG or long
#'  memory. More details on the precise simulation algorithm is available in the
#'  vignette.
#'@export
sim_UnivariateTrawl <-function(t, Delta=1, burnin=10,marginal =base::c("Poi", "NegBin"),trawl =base::c("Exp", "DExp", "supIG", "LM"),
                               v=0,m=0,theta=0, lambda1=0,lambda2=0,w=0,delta=0,gamma=0,alpha=0,H=0){

  T <- burnin + t #Simulate over [0,T] and later report the values on [burnin,T]

  #Generate the pairs {t_i, U_i}_{i=1,..,N_T}
  # First, determine the number of jumps in the interval [0,T]
  if(marginal=="NegBin"){ v<-m*base::abs(base::log(1-theta))}
  N_T = stats::rpois(1,v*T)

  #Second, simulate the jump times
  jumptimes = base::sort(stats::runif(N_T, min=0, max=T))


  #Third, simulate the jump heights of the Poisson basis
  jumpheights = stats::runif(N_T, min=0, max=1)

  #Fourth, draw the jump marks from the logarithmic series distribution
  if(marginal=="Poi"){Cs <-base::numeric(N_T)+1}
  if(marginal=="NegBin"){Cs <-Runuran::urlogarithmic(N_T,theta)}


  ####First, determine how many jumps happened up to each grid point and
  ####store the number of jumps in the following vector
  VectorOfJumpTimes <- base::vector(mode = "numeric", length = base::floor(T/Delta))
  cuttable <-base::table(base::cut(jumptimes,seq(0,T,Delta),include.lowest = TRUE))
  VectorOfJumpTimes[1]<- base::as.integer(cuttable[1])
  for(i in 2:base::floor(T/Delta)){
    VectorOfJumpTimes[i] <- VectorOfJumpTimes[i-1]+base::as.integer(cuttable[i])
  }



  ### Next, compute the trawl process
  TrawlProcess <- base::vector(mode = "numeric", length = base::floor(T/Delta))

  for(k in 1:base::floor(T/Delta))
  {
    NoJ <- VectorOfJumpTimes[k] #Number of jumps until time k*Delta
    if(NoJ>0){
      timediff <-k*Delta - jumptimes[1:NoJ]
      if(trawl=="Exp"){
        cond <-1-base::ceiling(jumpheights[1:NoJ]-trawl_Exp(-timediff,lambda1)) }#need jumpheights-function <0 to sum up
      if(trawl=="DExp"){
        cond <-1-base::ceiling(jumpheights[1:NoJ]-trawl_DExp(-timediff,w, lambda1, lambda2)) }
      if(trawl=="supIG"){
        cond <-1-base::ceiling(jumpheights[1:NoJ]-trawl_supIG(-timediff,delta,gamma)) }
      if(trawl=="LM"){
        cond <-1-base::ceiling(jumpheights[1:NoJ]-trawl_LM(-timediff,alpha,H)) }

      TrawlProcess[k] <- base::sum(cond*Cs[1:NoJ])
    }
  }

  #Cut off the  burn-in period
  bound1 = burnin/Delta
  bound2 = T/Delta
  length = bound2-bound1

  TP <- TrawlProcess[(bound1+1):bound2]

  return(TP)

}


#'Simulates a bivariate trawl process
#'@name sim_BivariateTrawl
#'@param t parameter which specifying the length of the time interval
#'  \eqn{[0,t]} for which a simulation of the trawl process is required.
#'@param Delta parameter \eqn{\Delta} specifying the length of the time step,
#'  the default is 1
#'@param burnin parameter specifying the length of the burn-in period at the
#'  beginning of the simulation
#'@param marginal parameter specifying the marginal distribution of the trawl
#'@param dependencetype Parameter specifying the type of dependence
#'@param trawl1 parameter specifying the type of the first trawl function
#'@param trawl2 parameter specifying the type of the second trawl function
#'@param v1,v2,v12 parameters of the Poisson distribution
#'@param kappa1,kappa2,kappa12,a1,a2 parameters of the (possibly bivariate)
#'  negative binomial distribution
#'@param lambda11,lambda12,w1 parameters of the exponential (or
#'  double-exponential) trawl function of the first process
#'@param delta1,gamma1 parameters of the supIG trawl function of the first
#'  process
#'@param alpha1,H1 parameter of the long memory trawl of the first process
#'@param lambda21,lambda22,w2 parameters of the exponential (or
#'  double-exponential) trawl function of the second process
#'@param delta2,gamma2 parameters of the supIG trawl function of the second
#'  process
#'@param alpha2,H2 parameter of the long memory trawl of the second process
#'@details This function simulates a bivariate trawl process with either Poisson
#'  or negative binomial marginal law. For the trawl function there are
#'  currently four choices: exponential, double-exponential, supIG or long
#'  memory. More details on the precise simulation algorithm is available in the
#'  vignette.
#'@export

sim_BivariateTrawl <-function(t, Delta=1, burnin=10,marginal =base::c("Poi", "NegBin"), dependencetype=base::c("fullydep","dep"),
                              trawl1 =base::c("Exp", "DExp", "supIG", "LM"),
                              trawl2 =base::c("Exp", "DExp", "supIG", "LM"),
                               v1=0,v2=0,v12=0,kappa1=0,kappa2=0,kappa12=0,a1=0,a2=0,
                              lambda11=0,lambda12=0,w1=0,delta1=0,gamma1=0,alpha1=0,H1=0,
                              lambda21=0,lambda22=0,w2=0,delta2=0,gamma2=0,alpha2=0,H2=0){
  T <- burnin + t #Simulate over [0,T] and later report the values on [burnin,T]

   if(dependencetype=="fullydep"){


    #Generate the pairs {t_i, U_i}_{i=1,..,N_T}
    # First, determine the number of jumps in the interval [0,T]
    if(marginal=="NegBin"){ v12<-kappa12*(base::log(1+a1+a2))}
    N_T = stats::rpois(1,v12*T)

    #Second, simulate the jump times
    jumptimes = base::sort(stats::runif(N_T, min=0, max=T))


    #Third, simulate the jump heights of the Poisson basis
    jumpheights = stats::runif(N_T, min=0, max=1)

    #Fourth, draw the jump marks from the logarithmic series distribution
    if(marginal=="Poi"){
      Cs1 <-base::numeric(N_T)+1
      Cs2 <-Cs1}
    if(marginal=="NegBin"){
      p1 <- a1/(a1+a2+1)
      p2 <- a2/(a1+a2+1)
      jumpmarks <- Bivariate_LSDsim(N_T,p1,p2)
      Cs1 <- jumpmarks[,1]
      Cs2 <- jumpmarks[,2]
     }


    ####First, determine how many jumps happened up to each grid point and
    ####store the number of jumps in the following vector
    VectorOfJumpTimes <- base::vector(mode = "numeric", length = base::floor(T/Delta))
    cuttable <-base::table(base::cut(jumptimes,seq(0,T,Delta),include.lowest = TRUE))
    VectorOfJumpTimes[1]<- base::as.integer(cuttable[1])
    for(i in 2:base::floor(T/Delta)){
      VectorOfJumpTimes[i] <- VectorOfJumpTimes[i-1]+base::as.integer(cuttable[i])
    }



  ### Next, compute the trawl process
  TrawlProcess1 <- base::vector(mode = "numeric", length = base::floor(T/Delta))
  TrawlProcess2 <- base::vector(mode = "numeric", length = base::floor(T/Delta))

  for(k in 1:base::floor(T/Delta))
    {
      NoJ <- VectorOfJumpTimes[k] #Number of jumps until time k*Delta
      if(NoJ>0){
        timediff <-k*Delta - jumptimes[1:NoJ]
        if(trawl1=="Exp"){
          cond1 <-1-base::ceiling(jumpheights[1:NoJ]-trawl_Exp(-timediff,lambda11)) }#need jumpheights-function <0 to sum up
        if(trawl1=="DExp"){
          cond1 <-1-base::ceiling(jumpheights[1:NoJ]-trawl_DExp(-timediff,w1, lambda11, lambda12)) }
        if(trawl1=="supIG"){
          cond1 <-1-base::ceiling(jumpheights[1:NoJ]-trawl_supIG(-timediff,delta1,gamma1)) }
        if(trawl1=="LM"){
          cond1 <-1-base::ceiling(jumpheights[1:NoJ]-trawl_LM(-timediff,alpha1,H1)) }

        if(trawl2=="Exp"){ #UPDATE PARAMETER NAMES FOR TWO PROCESSES
          cond2 <-1-base::ceiling(jumpheights[1:NoJ]-trawl_Exp(-timediff,lambda21)) }#need jumpheights-function <0 to sum up
        if(trawl2=="DExp"){
          cond2 <-1-base::ceiling(jumpheights[1:NoJ]-trawl_DExp(-timediff,w2, lambda21, lambda22)) }
        if(trawl2=="supIG"){
          cond2 <-1-base::ceiling(jumpheights[1:NoJ]-trawl_supIG(-timediff,delta2,gamma2)) }
        if(trawl2=="LM"){
          cond2 <-1-base::ceiling(jumpheights[1:NoJ]-trawl_LM(-timediff,alpha2,H2)) }

        TrawlProcess1[k] <- base::sum(cond1*Cs1[1:NoJ])
        TrawlProcess2[k] <- base::sum(cond2*Cs2[1:NoJ])
      }
    }
  #Cut off the  burn-in period
  bound1 = burnin/Delta
  bound2 = T/Delta
  length = bound2-bound1

  TP1 <- TrawlProcess1[(bound1+1):bound2]
  TP2 <- TrawlProcess2[(bound1+1):bound2]
  }
  if(dependencetype=="dep"){

    #Generate the pairs {t_i, U_i}_{i=1,..,N_T}
    # First, determine the number of jumps in the interval [0,T]
    if(marginal=="NegBin"){
      v12<-kappa12*(base::log(1+a1+a2))
      v1 <-kappa1*base::log(1+a1)
      v2 <-kappa2*base::log(1+a2)
    }
    N12_T <- stats::rpois(1,v12*T) #joint process
    if(N12_T==0){N12_T <- stats::rpois(1,v12*T)}

    if(v1==0){N1_T=0} else {N1_T <-stats::rpois(1,v1*T)}
    if(v2==0){N2_T=0} else {N2_T <-stats::rpois(1,v2*T)}

    #Second, simulate the jump times
    jumptimes12 <- base::sort(stats::runif(N12_T, min=0, max=T))
    jumptimes1 <- base::sort(stats::runif(N1_T, min=0, max=T))
    jumptimes2 <- base::sort(stats::runif(N2_T, min=0, max=T))

    alljumptimes1 <- base::sort(base::c(jumptimes12,jumptimes1))
    alljumptimes2 <- base::sort(base::c(jumptimes12,jumptimes2))


    #Third, simulate the jump heights of the Poisson basis
    jumpheights12 <- stats::runif(N12_T, min=0, max=1)
    jumpheights1  <- stats::runif(N1_T, min=0, max=1)
    jumpheights2  <- stats::runif(N2_T, min=0, max=1)

    #Fourth, draw the jump marks from the logarithmic series distribution
    if(marginal=="Poi"){
      Cs1 <-base::numeric(N1_T)+1
      Cs2 <-base::numeric(N2_T)+1
      Cs12_1 <-base::numeric(N12_T)+1
      Cs12_2 <-base::numeric(N12_T)+1
      }
    if(marginal=="NegBin"){
      p1 <- a1/(a1+a2+1)
      p2 <- a2/(a1+a2+1)
      jumpmarks12 <- Bivariate_LSDsim(N12_T,p1,p2)
      Cs12_1 <- jumpmarks12[,1]
      Cs12_2 <- jumpmarks12[,2]
      if(N1_T>0){Cs1 <- Runuran::urlogarithmic(N1_T,a1/(1+a1))} else{Cs1 <-0}
      if(N2_T>0){Cs2 <- Runuran::urlogarithmic(N2_T,a2/(1+a2))} else{Cs2 <-0}

    }


    ####First, determine how many jumps happened up to each grid point and
    ####store the number of jumps in the following vector
    VectorOfJumpTimes12 <- base::vector(mode = "numeric", length = floor(T/Delta))
    cuttable12 <-base::table(base::cut(jumptimes12,seq(0,T,Delta),include.lowest = TRUE))
    VectorOfJumpTimes12[1]<- base::as.integer(cuttable12[1])
    for(i in 2:base::floor(T/Delta)){
      VectorOfJumpTimes12[i] <- VectorOfJumpTimes12[i-1]+as.integer(cuttable12[i])
    }

    VectorOfJumpTimes1 <- base::vector(mode = "numeric", length = base::floor(T/Delta))
    cuttable1 <-base::table(base::cut(jumptimes1,seq(0,T,Delta),include.lowest = TRUE))
    VectorOfJumpTimes1[1]<- base::as.integer(cuttable1[1])
    for(i in 2:base::floor(T/Delta)){
      VectorOfJumpTimes1[i] <- VectorOfJumpTimes1[i-1]+base::as.integer(cuttable1[i])
    }

    VectorOfJumpTimes2 <- base::vector(mode = "numeric", length = base::floor(T/Delta))
    cuttable2 <-base::table(base::cut(jumptimes2,seq(0,T,Delta),include.lowest = TRUE))
    VectorOfJumpTimes2[1]<- base::as.integer(cuttable2[1])
    for(i in 2:base::floor(T/Delta)){
      VectorOfJumpTimes2[i] <- VectorOfJumpTimes2[i-1]+base::as.integer(cuttable2[i])
    }



    ### Next, compute the trawl process
    TrawlProcess1 <- base::vector(mode = "numeric", length = base::floor(T/Delta))
    TrawlProcess2 <- base::vector(mode = "numeric", length = base::floor(T/Delta))
    TrawlProcess12_1 <- base::vector(mode = "numeric", length = base::floor(T/Delta))
    TrawlProcess12_2 <- base::vector(mode = "numeric", length = base::floor(T/Delta))

    for(k in 1:base::floor(T/Delta))
    {
      NoJ <- VectorOfJumpTimes12[k] #Number of jumps until time k*Delta
      if(NoJ>0){
        timediff <-k*Delta - jumptimes12[1:NoJ]
        if(trawl1=="Exp"){ #UPDATE PARAMETER NAMES FOR TWO PROCESSES
          cond1 <-1-base::ceiling(jumpheights12[1:NoJ]-trawl_Exp(-timediff,lambda11)) }#need jumpheights-function <0 to sum up
        if(trawl1=="DExp"){
          cond1 <-1-base::ceiling(jumpheights12[1:NoJ]-trawl_DExp(-timediff,w1, lambda11, lambda12)) }
        if(trawl1=="supIG"){
          cond1 <-1-base::ceiling(jumpheights12[1:NoJ]-trawl_supIG(-timediff,delta1,gamma1)) }
        if(trawl1=="LM"){
          cond1 <-1-base::ceiling(jumpheights12[1:NoJ]-trawl_LM(-timediff,alpha1,H1)) }

        if(trawl2=="Exp"){ #UPDATE PARAMETER NAMES FOR TWO PROCESSES
          cond2 <-1-base::ceiling(jumpheights12[1:NoJ]-trawl_Exp(-timediff,lambda21)) }#need jumpheights-function <0 to sum up
        if(trawl2=="DExp"){
          cond2 <-1-base::ceiling(jumpheights12[1:NoJ]-trawl_DExp(-timediff,w2, lambda21, lambda22)) }
        if(trawl2=="supIG"){
          cond2 <-1-base::ceiling(jumpheights12[1:NoJ]-trawl_supIG(-timediff,delta2,gamma2)) }
        if(trawl2=="LM"){
          cond2 <-1-base::ceiling(jumpheights12[1:NoJ]-trawl_LM(-timediff,alpha2,H2)) }

        TrawlProcess12_1[k] <- base::sum(cond1*Cs12_1[1:NoJ])
        TrawlProcess12_2[k] <- base::sum(cond2*Cs12_2[1:NoJ])
      }

      NoJ <- VectorOfJumpTimes1[k] #Number of jumps until time k*Delta
      if(NoJ>0){
        timediff <-k*Delta - jumptimes1[1:NoJ]
        if(trawl1=="Exp"){ #UPDATE PARAMETER NAMES FOR TWO PROCESSES
          cond1 <-1-base::ceiling(jumpheights1[1:NoJ]-trawl_Exp(-timediff,lambda11)) }#need jumpheights-function <0 to sum up
        if(trawl1=="DExp"){
          cond1 <-1-base::ceiling(jumpheights1[1:NoJ]-trawl_DExp(-timediff,w1, lambda11, lambda12)) }
        if(trawl1=="supIG"){
          cond1 <-1-base::ceiling(jumpheights1[1:NoJ]-trawl_supIG(-timediff,delta1,gamma1)) }
        if(trawl1=="LM"){
          cond1 <-1-base::ceiling(jumpheights1[1:NoJ]-trawl_LM(-timediff,alpha1,H1)) }
        TrawlProcess1[k] <- base::sum(cond1*Cs1[1:NoJ])
      }
      NoJ <- VectorOfJumpTimes2[k] #Number of jumps until time k*Delta
      if(NoJ>0){
        timediff <-k*Delta - jumptimes2[1:NoJ]
        if(trawl2=="Exp"){
          cond2 <-1-base::ceiling(jumpheights2[1:NoJ]-trawl_Exp(-timediff,lambda21)) }#need jumpheights-function <0 to sum up
        if(trawl2=="DExp"){
          cond2 <-1-base::ceiling(jumpheights2[1:NoJ]-trawl_DExp(-timediff,w2, lambda21, lambda22)) }
        if(trawl2=="supIG"){
          cond2 <-1-base::ceiling(jumpheights2[1:NoJ]-trawl_supIG(-timediff,delta2,gamma2)) }
        if(trawl2=="LM"){
          cond2 <-1-base::ceiling(jumpheights2[1:NoJ]-trawl_LM(-timediff,alpha2,H2)) }
        TrawlProcess2[k] <- base::sum(cond2*Cs2[1:NoJ])
      }


    }
    #Cut off the  burn-in period
    bound1 = burnin/Delta
    bound2 = T/Delta
    length = bound2-bound1

    TP1 <- TrawlProcess1[(bound1+1):bound2]+TrawlProcess12_1[(bound1+1):bound2]
    TP2 <- TrawlProcess2[(bound1+1):bound2]+TrawlProcess12_2[(bound1+1):bound2]
  }







  TP <-base::cbind(TP1,TP2)

  return(TP)

}

#'Plots the bivariate histogram of two time series together with the univariate
#'histograms
#'@name plot_2and1hist
#'@param x vector of equidistant time series data
#'@param y vector of equidistant time series data (of the same length as x)
#'@details This function plots the bivariate histogram of two time series
#'  together with the univariate histograms
#'@return no return value
#'@export
plot_2and1hist <- function(x,y){
  if(base::NROW(x)==base::NROW(y)){
    df <- base::data.frame(x,y)
    h1 <- graphics::hist(df$x, plot=F)
    h2 <- graphics::hist(df$y,plot=F)
    top <- base::max(h1$counts, h2$counts)


    # Set the margins
    oldpar <- graphics::par()
    graphics::par(mar=base::c(3,3,1,1))
    graphics::layout(matrix(base::c(2,0,1,3),2,2,byrow=T),base::c(3,1), base::c(1,3))
    squash::hist2(df, base = 2, colFn=squash::rainbow2,key = squash::vkey, nz = 5, bty = 'l')

    graphics::par(mar=base::c(0,2,1,0))
    graphics::barplot(h1$counts, axes=T,  space=0, col='red')
    graphics::par(mar=base::c(2,0,0.5,1))
    graphics::barplot(h2$counts, axes=T, space=0, col='red', horiz=T)

    ###Back to default graphics setting
    graphics::par(mfrow=base::c(1,1))
    graphics::par(mar=base::c(5.1, 4.1, 4.1, 2.1), mgp=base::c(3, 1, 0), las=0)
  } else {
    print("Error: x and y do not have the same length")
  }

}


#'Plots the bivariate histogram of two time series together with the univariate
#'histograms using ggplot2
#'@name plot_2and1hist_gg
#'@param x vector of equidistant time series data
#'@param y vector of equidistant time series data (of the same length as x)
#'@param bivbins number of bins in the bivariate histogram
#'@param xbins number of bins in the histogram of x
#'@param ybins number of bins in the histogram of y
#'@details This function plots the bivariate histogram of two time series
#'  together with the univariate histograms
#'@return no return value
#'@export
plot_2and1hist_gg <- function(x,y, bivbins =50, xbins=30, ybins=30){
  if(base::NROW(x)==base::NROW(y)){
    df <- base::data.frame(x,y)
    g_biv <-ggplot2::ggplot(df, ggplot2::aes(x=x, y=y) ) +  ggplot2::geom_bin2d(bins = bivbins) +  ggplot2::scale_fill_continuous(type = "viridis") +  ggplot2::theme_bw()
    g_x <- ggplot2::qplot(x, geom="histogram")+ggplot2::xlab("x")+ggplot2::stat_bin(bins = xbins)
    g_y <- ggplot2::qplot(y, geom="histogram")+ggplot2::xlab("y")+ggplot2::stat_bin(bins = ybins)
    ggpubr::ggarrange(g_biv, ggpubr::ggarrange(g_x, g_y,  ncol = 2), nrow = 2, heights = c(2, 0.7))

  } else {
    print("Error: x and y do not have the same length")
  }

}
