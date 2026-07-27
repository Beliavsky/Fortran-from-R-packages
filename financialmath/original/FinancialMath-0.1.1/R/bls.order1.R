bls.order1 = function(S,K,r,t,sd,D=0){

  ##Error Handling
  test=c(S,K,r,t,sd)
  if(t < 0)stop('Error: square root of negative number \n  varibale "t" must be positive')
  if(any(test==0))stop('Error: divides by 0 \n  only dividend yield (D) can be 0')
  if(sd<0)warning('volitility (sd) input as negative')
  if(is.numeric(test)==F) stop("inputs must be numeric")
  if((S<0 & K>0)|((S>0 & K<0)))stop("Error: log of negative number \n  S & K must both be positive or negative")

  ##d1 & d2
  d1 = (log(S/K)+(r-D+.5*sd^2)*t)/(sd*sqrt(t))
  d2 = (log(S/K)+(r-D-.5*sd^2)*t)/(sd*sqrt(t))
  ##Greek Values
  Call_Value = S*exp(-D*t)*pnorm(d1)-K*exp(-r*t)*pnorm(d2)
  Put_Value = -S*exp(-D*t)*pnorm(-d1)+K*exp(-r*t)*pnorm(-d2)
  Call_Delta = exp(-D*t)*pnorm(d1)
  Put_Delta = exp(-D*t)*(pnorm(d1)-1)
  Call_Theta = (-sd*S*exp(-D*t)*dnorm(d1))/(2*sqrt(t))+D*S*pnorm(-d1)*exp(-D*t)-
    r*K*exp(-r*t)*pnorm(d2)
  Put_Theta = (-sd*S*exp(-D*t)*dnorm(-d1))/(2*sqrt(t))-D*S*pnorm(d1)*exp(-D*t)+
    r*K*exp(-r*t)*pnorm(-d2)
  Vega = S*sqrt(t)*exp(-D*t)*dnorm(d1)

  ##Output
  return(matrix(c(Call_Value,Put_Value,Call_Delta,Put_Delta,Call_Theta,Put_Theta,Vega,Vega),nrow=4,
                ncol=2, byrow=T,
                dimnames=list(c("Price","Delta","Theta","Vega"),
                              c("Call","Put"))))
}

