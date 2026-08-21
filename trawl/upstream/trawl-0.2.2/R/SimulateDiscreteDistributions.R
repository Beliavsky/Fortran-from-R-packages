#'Simulates from the bivariate negative binomial distribution
#'@name Bivariate_NBsim
#'@param N number of data points to be simulated
#'@param kappa parameter \eqn{\kappa} of the bivariate negative binomial distribution
#'@param p1 parameter \eqn{p_1} of the bivariate negative binomial distribution
#'@param p2 parameter \eqn{p_2} of the bivariate negative binomial distribution
#'@return An \eqn{N\times 2} matrix with \eqn{N} simulated values from the bivariate negative
#'  binomial distribution
#'@details A random vector \eqn{{\bf X}=(X_1,X_2)'} is said to follow the
#'  bivariate negative binomial distribution with parameters \eqn{\kappa, p_1,
#'  p_2} if its probability mass function is given by \deqn{ P({\bf X}={\bf
#'  x})=\frac{\Gamma(x_1+x_2+\kappa)}{x_1!x_2!
#'  \Gamma(\kappa)}p_1^{x_1}p_2^{x_2}(1-p_1-p_2)^{\kappa},} where,
#'  for \eqn{i=1,2}, \eqn{x_i\in\{0,1,\dots\}}, \eqn{0<p_i<1} such that
#'  \eqn{p_1+p_2<1} and \eqn{\kappa>0}.

#'@export

Bivariate_NBsim <-function(N, kappa, p1, p2){

  #First check that N is a positive integer and that kappa, p1 and p2 satisfy the parameter constraints
  if((N>0)&&(N%%1==0)&&(p1>0)&&(p2>0)&&(p1<1)&&(p2<1)&&((p1+p2)<1)&&(kappa>0)&&(kappa<Inf)){

    #Simulate X_1 from the univariate negative Binomial distribution
    X1 <- stats::rnbinom(N,size=kappa,prob=1-p1/(1-p2))

    X2 <- numeric(N)
    for(i in 1:N){
      #Simulate X_2|X_1 from the univariate negative Binomial distribution
      X2[i] <- stats::rnbinom(1,size=kappa+X1[i],prob=1-p2)
    }

    X <- base::cbind(X1,X2)
    X

  }

  else{    print("Invalid parameter choice")  }

}


#'Simulates from the bivariate logarithmic series distribution
#'@name Bivariate_LSDsim
#'@param N number of data points to be simulated
#'@param p1 parameter \eqn{p_1} of the bivariate logarithmic series distribution
#'@param p2 parameter \eqn{p_2} of the bivariate logarithmic series distribution
#'@return An \eqn{N \times 2} matrix with \eqn{N} simulated values from the
#'  bivariate logarithmic series distribution
#'@details The probability mass function of a random vector \eqn{X=(X_1,X_2)' }
#'  following the bivariate logarithmic series distribution with parameters
#'  \eqn{0<p_1, p_2<1} with \eqn{p:=p_1+p_2<1} is given by
#'  \deqn{P(X_1=x_1,X_2=x_2)=\frac{\Gamma(x_1+x_2)}{x_1!x_2!}
#'  \frac{p_1^{x_1}p_2^{x_2}}{(-\log(1-p))},} for \eqn{x_1,x_2=0,1,2,\dots} such
#'  that \eqn{x_1+x_2>0}. The simulation proceeds in two steps: First, \eqn{X_1}
#'  is simulated from the modified logarithmic distribution with parameters
#'  \eqn{\tilde p_1=p_1/(1-p_2)} and \eqn{\delta_1=\log(1-p_2)/\log(1-p)}. Then
#'  we simulate \eqn{X_2} conditional on \eqn{X_1}. We note that
#'  \eqn{X_2|X_1=x_1} follows the  logarithmic series distribution with
#'  parameter \eqn{p_2} when \eqn{x_1=0}, and the  negative binomial
#'  distribution with parameters \eqn{(x_1,p_2)} when  \eqn{x_1>0}.
#'@export
Bivariate_LSDsim <-function(N, p1, p2){

  #First check that N is a positive integer and that p1 and p2 satisfy the parameter constraints
  if((N>0)&&(N%%1==0)&&(p1>0)&&(p2>0)&&(p1<1)&&(p2<1)&&((p1+p2)<1)){

    #Compute the parameter values of the modified LSD
    delta1 <- base::log(1-p2)/log(1-p1-p2)
    tilde_p1 <- p1/(1-p2)

    #Simulate X1 from the modified logarithmic distribution
    Log <- Runuran::urlogarithmic(N,tilde_p1) #rlog(N,tilde_p1)
    Bernoulli <- stats::rbinom(N,1,1-delta1)
    X1 <- Log * Bernoulli
    X2 <- numeric(N)

    for(i in 1:N){
      x1 <- X1[i]
      if(x1==0){ X2[i] <- Runuran::urlogarithmic(1,p2)}# rlog(1,p2) }
      if(x1>0){  X2[i] <- stats::rnbinom(1,size=x1,prob=1-p2) }
    }

    X <- cbind(X1,X2)
    X

  }

  else{    print("Invalid parameter choice")  }

}


#'Simulates from the trivariate logarithmic series distribution
#'@name Trivariate_LSDsim
#'@param N number of data points to be simulated
#'@param p1 parameter \eqn{p1} of the trivariate logarithmic series distribution
#'@param p2 parameter \eqn{p2} of the trivariate logarithmic series distribution
#'@param p3 parameter \eqn{p3} of the trivariate logarithmic series distribution
#'@return An \eqn{N \times 3} matrix with \eqn{N} simulated values from the
#'  trivariate logarithmic series distribution
#'@details The probability mass function of a random vector
#'  \eqn{X=(X_1,X_2,X_3)' } following the trivariate logarithmic series
#'  distribution with parameters \eqn{0<p_1, p_2, p_3<1} with
#'  \eqn{p:=p_1+p_2+p_3<1} is given by
#'  \deqn{P(X_1=x_1,X_2=x_2,X_3=x_3)=\frac{\Gamma(x_1+x_2+x_3)}{x_1!x_2!x_3!}
#'  \frac{p_1^{x_1}p_2^{x_2}p_3^{x_3}}{(-\log(1-p))},} for
#'  \eqn{x_1,x_2,x_3=0,1,2,\dots} such that \eqn{x_1+x_2+x_3>0}.
#'
#'  The simulation proceeds in two steps: First, \eqn{X_1} is simulated from the
#'  modified logarithmic distribution with parameters \eqn{\tilde
#'  p_1=p_1/(1-p_2-p_3)} and \eqn{\delta_1=\log(1-p_2-p_3)/\log(1-p)}. Then we
#'  simulate \eqn{(X_2,X_3)'} conditional on \eqn{X_1}. We note that
#'  \eqn{(X_2,X_3)'|X_1=x_1} follows the bivariate logarithmic series
#'  distribution with parameters \eqn{(p_2,p_3)} when \eqn{x_1=0}, and the
#'  bivariate negative binomial distribution with parameters \eqn{(x_1,p_2,p_3)}
#'  when  \eqn{x_1>0}.
#'@export


Trivariate_LSDsim <-function(N, p1, p2, p3){

  #First check that N is a positive integer and that p1 and p2 satisfy the parameter constraints
  if((N>0)&&(N%%1==0)&&(p1>0)&&(p2>0)&&(p3>0)&&(p1<1)&&(p2<1)&&(p3<1)&&((p1+p2+p3)<1)){

    #Compute the parameter values of the modified LSD
    delta1 <- base::log(1-p2-p3)/log(1-p1-p2-p3)
    tilde_p1 <- p1/(1-p2-p3)

    #Simulate X1 from the modified logarithmic distribution
    Log <- Runuran::urlogarithmic(N,tilde_p1) #rlog(N,tilde_p1)
    Bernoulli <- stats::rbinom(N,1,1-delta1)
    X1 <- Log * Bernoulli
    X2 <- numeric(N)
    X3 <- numeric(N)

    for(i in 1:N){
      x1 <- X1[i]
      if(x1==0){
        Y <- Bivariate_LSDsim(1,p2,p3)
        X2[i] <- Y[1]
        X3[i] <- Y[2]
      }
      if(x1>0){
        Y <- Bivariate_NBsim(1,x1,p2,p3)
        X2[i] <- Y[1]
        X3[i] <- Y[2]
      }
    }

    X <- cbind(X1,X2,X3)
    X

  }

  else{    print("Invalid parameter choice")  }

}






#'Computes the mean of the logarithmic series distribution
#'@name LSD_Mean
#'@param p parameter of the logarithmic series distribution
#'@return Mean of the logarithmic series distribution
#'@details A random variable \eqn{X} has logarithmic series distribution with
#'  parameter \eqn{0<p<1}  if \deqn{P(X=x)=(-1)/(\log(1-p))p^x/x, \mbox{ for }
#'  x=1,2,\dots.}
#'@export
LSD_Mean <-function(p){p/((1-p)*(-log(1-p)))}


#'Computes the variance of the logarithmic series distribution
#'@name LSD_Var
#'@param p parameter of the logarithmic series distribution
#'@return Variance of the logarithmic series distribution
#'@details A random variable \eqn{X} has logarithmic series distribution with
#'  parameter \eqn{0<p<1}  if \deqn{P(X=x)=(-1)/(\log(1-p))p^x/x, \mbox{ for }
#'  x=1,2,\dots.}
#'@export

LSD_Var<-function(p){(-p)*(p+log(1-p))/((1-p)*log(1-p))^2}


#'Computes the mean of the modified logarithmic series distribution
#'@name ModLSD_Mean
#'@param delta parameter \eqn{\delta} of the modified logarithmic series
#'  distribution
#'@param p parameter of the modified logarithmic series distribution
#'@return Mean of the modified logarithmic series distribution
#'@details A random variable \eqn{X} has modified logarithmic series
#'  distribution with parameters \eqn{0 \le \delta <1} and \eqn{0<p<1}  if
#'  \eqn{P(X=0)=(1-\delta)} and \deqn{P(X=x)=(1-\delta)(-1)/(\log(1-p))p^x/x,
#'  \mbox{ for  } x=1,2,\dots.}
#'@export
ModLSD_Mean<-function(delta,p){(1-delta)*LSD_Mean(p)}


#'Computes the variance of the modified logarithmic series distribution
#'@name ModLSD_Var
#'@param delta parameter \eqn{\delta} of the modified logarithmic series
#'  distribution
#'@param p parameter of the modified logarithmic series distribution
#'@return Mean of the modified logarithmic series distribution
#'@details A random variable \eqn{X} has modified logarithmic series
#'  distribution with parameters \eqn{0\le \delta <1} and \eqn{0<p<1}  if
#'  \eqn{P(X=0)=(1-\delta)} and \deqn{P(X=x)=(1-\delta)(-1)/(\log(1-p))p^x/x,
#'  \mbox{ for }  x=1,2,\dots.}
#'@export
ModLSD_Var<-function(delta,p){
  EB <- (1-delta)
  EB2 <-delta*(1-delta)+(1-delta)^2
  VarB<-delta*(1-delta)
  VarAgivenB <-LSD_Var(p)
  EAgivenB <- LSD_Mean(p)
  R <- EB2*VarAgivenB+VarB*EAgivenB^2
  return(R)

}


#'Computes the correlation of the components of a bivariate vector following the
#'bivariate  logarithmic series distribution
#'@name BivLSD_Cor
#'@param p1 parameter \eqn{p_1} of the bivariate  logarithmic series
#'  distribution
#'@param p2 parameter \eqn{p_2} of the bivariate  logarithmic series
#'  distribution
#'@return Correlation of the components of a bivariate vector following the
#'  bivariate  logarithmic series distribution
#'@export

BivLSD_Cor <-function(p1,p2){
  gamma <- 1-p1-p2
  L<- -base::log(gamma)
  correlation <- base::sign(L-1)*base::sqrt(p1*p2*(L-1)^2/((p1*(L-1)+gamma*L)*(p2*(L-1)+gamma*L)))
  return(correlation )
}


#'Computes the covariance of the components of a bivariate vector following the
#'bivariate  logarithmic series distribution
#'@name BivLSD_Cov
#'@param p1 parameter \eqn{p_1} of the bivariate  logarithmic series
#'  distribution
#'@param p2 parameter \eqn{p_2} of the bivariate  logarithmic series
#'  distribution
#'@return Covariance of the components of a bivariate vector following the
#'  bivariate  logarithmic series distribution
#'@export

BivLSD_Cov <-function(p1,p2){

  tildep1 <-p1/(1-p2)
  tildep2 <-p2/(1-p1)


  delta1<-base::log(1-p2)/base::log(1-p1-p2)
  delta2<-base::log(1-p1)/base::log(1-p1-p2)

  return(BivLSD_Cor(p1,p2)*(base::sqrt(ModLSD_Var(delta1, tildep1)*ModLSD_Var(delta2, tildep2))))
}


#'Computes the covariance of the components of a bivariate vector following the
#'bivariate modified  logarithmic series distribution
#'@name BivModLSD_Cov
#'@param delta parameter \eqn{\delta} of the bivariate modified   logarithmic
#'  series distribution
#'@param p1 parameter \eqn{p_1} of the bivariate modified   logarithmic series
#'  distribution
#'@param p2 parameter \eqn{p_2} of the bivariate  modified  logarithmic series
#'  distribution
#'@return Covariance of the components of a bivariate vector following the
#'  bivariate  modified  logarithmic series distribution
#'@export

BivModLSD_Cov <-function(delta,p1,p2){
  EB <- (1-delta)
  EB2 <-delta*(1-delta)+(1-delta)^2
  VarB<-delta*(1-delta)

  tildep1 <-p1/(1-p2)
  tildep2 <-p2/(1-p1)


  delta1<-base::log(1-p2)/base::log(1-p1-p2)
  delta2<-base::log(1-p1)/base::log(1-p1-p2)

  EA1givenB <- (1-delta1)*LSD_Mean(tildep1)
  EA2givenB <- (1-delta2)*LSD_Mean(tildep2)


  covariance <-EB2*BivLSD_Cov(p1,p2)+EA1givenB*EA2givenB*VarB
  return(covariance)
}

#'Computes the correlation of the components of a bivariate vector following the
#'bivariate modified  logarithmic series distribution
#'@name BivModLSD_Cor
#'@param delta parameter \eqn{\delta} of the bivariate modified   logarithmic
#'  series distribution
#'@param p1 parameter \eqn{p_1} of the bivariate modified   logarithmic series
#'  distribution
#'@param p2 parameter \eqn{p_2} of the bivariate  modified  logarithmic series
#'  distribution
#'@return Covariance of the components of a bivariate vector following the
#'  bivariate  modified  logarithmic series distribution
#'@export

BivModLSD_Cor <-function(delta,p1,p2){
  tildep1 <-p1/(1-p2)
  tildep2 <-p2/(1-p1)


  delta1<-base::log(1-p2)/base::log(1-p1-p2)
  delta2<-base::log(1-p1)/base::log(1-p1-p2)

  covariance <-BivModLSD_Cov(delta,p1,p2)/(base::sqrt(ModLSD_Var(delta1, tildep1)*ModLSD_Var(delta2, tildep2)))
  return(covariance)
}
