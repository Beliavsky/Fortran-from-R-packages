#' This function constructs GJRSK model of given data and parameters.
#' @param params vector of GJRSK model parameters(p1,const2,p2,q2,r2,const3,p3,q3,r3,const4,p4,q4,r4)
#' @param data vector time series data
#'
#' @return list of conditional mean(mu), variance(h), skewness(sk) and kurtosis(ku)

gjrsk_construct<-function(params,data){
  T<-length(data)

  mu <- rep(base::mean(data),T)
  h <- rep(stats::var(data),T)
  sk <- rep(skewness(data),T)
  ku <- rep(kurtosis(data),T)

  para_m<-params[1]
  para_h<-params[2:5]
  para_sk<-params[6:9]
  para_ku<-params[10:13]

  for(t in 2:T){
    mu[t]<-para_m * data[t-1]
    ut<-(data[t-1]-mu[t-1])
    eta<-(data[t-1]-mu[t-1])/sqrt(h[t-1])
    #Indicator Function
    I<-ifelse(ut<0,1,0)
    h[t]<-para_h %*% c(1,ut^2,(ut^2)*I,h[t-1])
    sk[t]<-para_sk %*% c(1,eta^3,(eta^3)*I,sk[t-1])
    ku[t]<-para_ku %*% c(1,eta^4,(eta^4)*I,ku[t-1])
  }

  return(list(mu=mu,h=h,sk=sk,ku=ku))

}

#' This function calculates the log-likelihood of GJRSK model.
#' @param params vector of GJRSK model parameters(p1,const2,p2,q2,r2,const3,p3,q3,r3,const4,p4,q4,r4)
#' @param data vector time series data
#'
#' @return (negative) log-likelihood of GJRSK model
gjrsk_lik<-function(params,data){
  T<-length(data)

  t <- 2:T
  likelihoods <- numeric(0)

  GJRSK <- gjrsk_construct(params,data)
  mu = GJRSK$mu[t]
  h = GJRSK$h[t]
  sk = GJRSK$sk[t]
  ku = GJRSK$ku[t]

  eta<-(data[t]-mu)/sqrt(h)

  f <- log((1 + (sk/6)*(eta^3 - 3*eta) + ((ku-3)/24)*(eta^4 - 6*eta^2+3))^2)
  g <- log(1 + (sk^2)/6 + ((ku-3)^2)/24)

  likelihoods <- -0.5*(log(h) + eta^2 ) + f - g
  #likelihoods <- -0.5*(log(h) + eta^2 + log(2*pi)) + f - g
  likelihoods <- - likelihoods

  LLF <-sum(likelihoods)

  if( is.nan(LLF) ){
    LLF <- 1e+06
  }

  return(LLF)
}

#' This function is inequality equation of GJRSK parameters used in optimization process(Rsolnp).
#' @param params vector of GJRSK model parameters(p1,const2,p2,q2,r2,const3,p3,q3,r3,const4,p4,q4,r4)
#' @param data vector time series data
#'
#' @return  upper bound >parameters > lower bound
gjrsk_ineqfun<-function(params,data){
  # -1 < a1 <1,
  # b0>0, b1>0, b3>0, and 0< b1 + b2 + b3<1,
  # -1 < c1 + c2 + c3< 1,
  # d0>0, d1>0, d3>0, and 0< d1 + d2 + d3<1

  para_m<-params[1]
  para_h<-params[2:5]
  para_sk<-params[6:9]
  para_ku<-params[10:13]

  return(c(para_m,para_h,para_sk[-1],para_ku,sum(para_h[-1]),sum(para_sk[-1]),sum(para_ku[-1])))

}


#' This function estimates GJRSK model's parameters.
#' @param data vector time series data
#'
#' @return list of parameters,standard errors of parameters,t-statistics,the minimum value of log-likelihood,AIC and BIC.
gjrsk_est<-function(data){
  X  <- data[-length(data)]
  Y  <- data[-1]
  a1 <- stats::cov(Y,X)/stats::var(X)
  b1 <- 0.01
  b2 <- 0.1
  b3 <- 0.85
  b0 <- (1-(sum(b1)+sum(b2)))*stats::var(data)
  c1 <- 0.01
  c2 <- 0.1
  c3 <- 0.7
  c0 <- (1-(sum(c1)+sum(c2)))*skewness(data)
  d1 <- 0.01
  d2 <- 0.05
  d3 <- 0.80
  d0 <- (1-(sum(d1)+sum(d2)))*(kurtosis(data))
  init <-c(a1,b0,b1,b2,b3,c0,c1,c2,c3,d0,d1,d2,d3)

  # -1 < a1 <1,
  # b0>0, b1>0, b3>0, and 0< b1 + b2 + b3<1,
  # -1 < c1 + c2 + c3< 1,
  # d0>0, d1>0, d3>0, and 0< d1 + d2 + d3<1

  aLB<-c(-1)
  bLB<-rep(0,4)
  cLB<-rep(-1,3)
  dLB<-rep(0,4)
  sumbLB<-c(0)
  sumcLB<-c(-1)
  sumdLB<-c(0)
  ineqLB<-c(aLB,bLB,cLB,dLB,sumbLB,sumcLB,sumdLB)

  aUB<-c(1)
  bUB<-c(Inf,rep(1,3))
  cUB<-rep(1,3)
  dUB<-c(Inf,rep(1,3))
  sumbUB<-c(1)
  sumcUB<-c(1)
  sumdUB<-c(1)
  ineqUB<-c(aUB,bUB,cUB,dUB,sumbUB,sumcUB,sumdUB)

  sol<-Rsolnp::solnp(pars=init , fun=gjrsk_lik, ineqfun = gjrsk_ineqfun, ineqLB = ineqLB, ineqUB = ineqUB, data=data)
  #sol<-solnp(pars=init , fun=gjrsk_lik, ineqfun = gjrsk_ineqfun, ineqLB = ineqLB, ineqUB = ineqUB,control=list(tol=1e-3), data=data)
  params<-sol$par
  loglik<-sol$values
  loglik<-loglik[length(loglik)]

  hessian<-sol$hessian[1:length(params),1:length(params)]
  #hessian(func=gjrsk_lik, x=params, data=data)
  stderrors <- sqrt(1/length(data)*diag(hessian))
  tstats <- params/stderrors
  AIC <- -2*loglik + 2*length(params)
  BIC <- -2*loglik + length(params)*log(length(data))

  return(list(params=params,stderrors=stderrors,tstats=tstats,loglik=loglik,AIC=AIC,BIC=BIC))

}


#' This function forcasts conditional mean,variance,skewness and kurtosis with given GJRSK model.
#' @param params vector of GJRSK model parameters(p1,const2,p2,q2,r2,const3,p3,q3,r3,const4,p4,q4,r4)
#' @param data vector time series data
#' @param max_forecast how long does this function forecast(Default value is 20)
#'
#' @return list of predicted conditional mean,variance,skewness and kurtosis
gjrsk_fcst<-function(params,data,max_forecast=20){
  T<-length(data)

  t <- 2:T

  GJRSK <- gjrsk_construct(params,data)
  mu = GJRSK$mu[t]
  h = GJRSK$h[t]
  sk = GJRSK$sk[t]
  ku = GJRSK$ku[t]

  para_m<-params[1]
  para_h<-params[2:5]
  para_sk<-params[6:9]
  para_ku<-params[10:13]

  T<-length(h)
  mu_fcst<-para_m * data[T]
  ut<-(data[T]-mu[T])
  eta<-(data[T]-mu[T])/sqrt(h[T])
  #Indicator Function
  I<-ifelse(ut<0,1,0)
  h_fcst<-para_h %*% c(1,ut^2,(ut^2)*I,h[T])
  sk_fcst<-para_sk %*% c(1,eta^3,(eta^3)*I,sk[T])
  ku_fcst<-para_ku %*% c(1,eta^4,(eta^4)*I,ku[T])


  for(i in 2:max_forecast){
    mu_fcst<-c(h_fcst, para_m * h_fcst[i-1])
    h_fcst<-c(h_fcst, para_h[1] + sum(para_h[-1] * h_fcst[i-1]))
    sk_fcst<-c(sk_fcst, para_sk[1] + sum(para_sk[-1] * sk_fcst[i-1]))
    ku_fcst<-c(ku_fcst, para_ku[1] + sum(para_ku[-1] * ku_fcst[i-1]))

  }

  return(list(mu_fcst=mu_fcst,h_fcst=h_fcst,sk_fcst=sk_fcst,ku_fcst=ku_fcst))

}
