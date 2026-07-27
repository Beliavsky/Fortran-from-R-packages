#' This function constructs GARCHSK model of given data and parameters.
#' @param params vector of GJRSK model parameters(p1,const2,p2,q2,const3,p3,q3,const4,p4,q4)
#' @param data vector time series data
#'
#' @return list of conditional mean(mu), variance(h), skewness(sk) and kurtosis(ku)

garchsk_construct<-function(params,data){
  T<-length(data)

  mu <- rep(base::mean(data),T)
  h <- rep(stats::var(data),T)
  sk <- rep(skewness(data),T)
  ku <- rep(kurtosis(data),T)

  para_m<-params[1]
  para_h<-params[2:4]
  para_sk<-params[5:7]
  para_ku<-params[8:10]

  for(t in 2:T){
    mu[t]<-para_m * data[t-1]
    h[t]<-para_h %*% c(1,(data[t-1]-mu[t-1])^2,h[t-1])
    sk[t]<-para_sk %*% c(1,((data[t-1]-mu[t-1])/sqrt(h[t-1]))^3,sk[t-1])
    ku[t]<-para_ku %*% c(1,((data[t-1]-mu[t-1])/sqrt(h[t-1]))^4,ku[t-1])
  }

  return(list(mu=mu,h=h,sk=sk,ku=ku))

}

#' This function calculates the log-likelihood of GARCHSK model.
#' @param params vector of GARCHSK model parameters(p1,const2,p2,q2,const3,p3,q3,const4,p4,q4)
#' @param data vector time series data
#'
#' @return (negative) log-likelihood of GJRSK model
garchsk_lik<-function(params,data){
  T<-length(data)

  t <- 2:T
  likelihoods <- numeric(0)

  GARCHSK <- garchsk_construct(params,data)
  mu = GARCHSK$mu[t]
  h = GARCHSK$h[t]
  sk = GARCHSK$sk[t]
  ku = GARCHSK$ku[t]

  std<-(data[t]-mu)/sqrt(h)

  f <- log((1 + (sk/6)*(std^3 - 3*std) + ((ku-3)/24)*(std^4 - 6*(std^2)+3))^2)
  g <- log(1 + (sk^2)/6 + ((ku-3)^2)/24)

  likelihoods <- -0.5*(log(h) + std^2 ) + f - g
  #likelihoods <- -0.5*(log(h) + std^2 + log(2*pi)) + f - g
  likelihoods <- - likelihoods

  LLF <-sum(likelihoods)

  if( is.nan(LLF) ){
    LLF <- 1e+06
  }

  return(LLF)
}


#' This function is inequality equation of GARCHSK parameters used in optimization process(Rsolnp).
#' @param params vector of GARCHSK model parameters(p1,const2,p2,q2,r2,const3,p3,q3,r3,const4,p4,q4,r4)
#' @param data vector time series data
#'
#' @return  upper bound >parameters > lower bound
garchsk_ineqfun<-function(params,data){
  # -1 < a1 <1,
  # b0>0, b1>0, b2>0, and 0< b1 + b2 <1,
  # -1 < c1 + c2 < 1, d0>0, d1>0, d2>0, and 0< d1+d2< 1


  para_m<-params[1]
  para_h<-params[2:4]
  para_sk<-params[5:7]
  para_ku<-params[8:10]
  return(c(para_m,para_h,para_sk[-1],para_ku,sum(para_h[-1]),sum(para_sk[-1]),sum(para_ku[-1])))

}


#' This function estimates GARCHSK model's parameters.
#' @param data vector time series data
#'
#' @return list of parameters,standard errors of parameters,t-statistics,the minimum value of log-likelihood,AIC and BIC.
garchsk_est<-function(data){
  X  <- data[-length(data)]
  Y  <- data[-1]
  a1 <- stats::cov(Y,X)/stats::var(X)
  b1 <- 0.01
  b2 <- 0.90
  b0 <- (1-(b1+b2))*stats::var(data)
  c1 <- 0.01
  c2 <- 0.7
  c0 <- (1-(c1+c2))*skewness(data)
  d1 <- 0.01
  d2 <- 0.80
  d0 <- (1-(d1+d2))*(kurtosis(data))
  init <-c(a1,b0,b1,b2,c0,c1,c2,d0,d1,d2)

  # -1 < a1 <1,
  # b0>0, b1>0, b2>0, and 0< b1 + b2 <1,
  # -1 < c1 + c2 < 1, d0>0, d1>0, d2>0, and 0< d1+d2< 1

  aLB<-c(-1)
  bLB<-c(0,rep(0,length(b1)),rep(0,length(b2)))
  cLB<-c(rep(-1,length(c1)),rep(-1,length(c2)))
  dLB<-c(0,rep(0,length(d1)),rep(0,length(d2)))
  sumbLB<-c(0)
  sumcLB<-c(-1)
  sumdLB<-c(0)
  ineqLB<-c(aLB,bLB,cLB,dLB,sumbLB,sumcLB,sumdLB)

  aUB<-c(1)
  bUB<-c(Inf,rep(1,length(b1)),rep(1,length(b2)))
  cUB<-c(rep(1,length(c1)),rep(1,length(c2)))
  dUB<-c(Inf,rep(1,length(d1)),rep(1,length(d2)))
  sumbUB<-c(1)
  sumcUB<-c(1)
  sumdUB<-c(1)
  ineqUB<-c(aUB,bUB,cUB,dUB,sumbUB,sumcUB,sumdUB)

  sol<-Rsolnp::solnp(pars=init , fun=garchsk_lik, ineqfun = garchsk_ineqfun, ineqLB = ineqLB, ineqUB = ineqUB, data=data)
  params<-sol$par
  loglik<-sol$values
  loglik<-loglik[length(loglik)]

  hessian<-sol$hessian[1:length(params),1:length(params)]
  #hessian(func=garchsk_lik, x=paraA, data=data,p1=p1,p2=p2,p3=p3,p4=p4,q2=q2,q3=q3,q4=q4)
  stderrors <- sqrt(1/length(data)*diag(hessian))
  tstats <- params/stderrors
  AIC <- -2*loglik + 2*length(params)
  BIC <- -2*loglik + length(params)*log(length(data))

  return(list(params=params,stderrors=stderrors,tstats=tstats,loglik=loglik,AIC=AIC,BIC=BIC))

}

#' This function forcasts conditional mean,variance,skewness and kurtosis with given GARCHSK model.
#' @param params vector of GARCHSK model parameters(p1,const2,p2,q2,const3,p3,q3,const4,p4,q4)
#' @param data vector time series data
#' @param max_forecast how long does this function forecast(Default value is 20)
#'
#' @return list of predicted conditional mean,variance,skewness and kurtosis
garchsk_fcst<-function(params,data,max_forecast=20){
  T<-length(data)
  t <- 2:T

  GARCHSK <- garchsk_construct(params,data)
  mu = GARCHSK$mu[t]
  h = GARCHSK$h[t]
  sk = GARCHSK$sk[t]
  ku = GARCHSK$ku[t]

  para_m<-params[1]
  para_h<-params[2:4]
  para_sk<-params[5:7]
  para_ku<-params[8:10]

  T<-length(h)
  mu_fcst<-para_m * data[T]
  h_fcst<-para_h %*% c(1,(data[T]-mu[T])^2,h[T])
  sk_fcst<-para_sk %*% c(1,((data[T]-mu[T])/sqrt(h[T]))^3,sk[T])
  ku_fcst<-para_ku %*% c(1,((data[T]-mu[T])/sqrt(h[T]))^4,ku[T])


  for(i in 2:max_forecast){
    mu_fcst<-c(h_fcst, para_m * h_fcst[i-1])
    h_fcst<-c(h_fcst, para_h[1] + sum(para_h[-1] * h_fcst[i-1]))
    sk_fcst<-c(sk_fcst, para_sk[1] + sum(para_sk[-1] * sk_fcst[i-1]))
    ku_fcst<-c(ku_fcst, para_ku[1] + sum(para_ku[-1] * ku_fcst[i-1]))

  }

  return(list(mu_fcst=mu_fcst,h_fcst=h_fcst,sk_fcst=sk_fcst,ku_fcst=ku_fcst))

}
