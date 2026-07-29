RM2006 = function(data, tau0 = 1560, tau1 = 4, kmax = 14, rho = sqrt(2)){
T = dim(data)[[1]]
K = dim(data)[[2]]
temp = array(0,dim = c(K,K,T))
Ht = Httilde = array(0,dim = c(K,K,T+1))


for (t in 1:T){
  temp[,,t] = tcrossprod(data[t,])
}  
data = temp

tauks = tau1*rho^((1:kmax)-1)
w = 1-log(tauks)/log(tau0)  
w = w/sum(w)

for (k in 1:kmax){
  mu = exp(-1/tauks[k])
  endPoint = max(min(floor(log(0.01)/log(mu)),T),k)
  weights = (1-mu)*mu^(0:(endPoint-1))
  weights = weights/sum(weights)
  backCast = matrix(0,K,K)
  for (i in 1:endPoint){
    backCast = backCast + weights[i]*data[,,i]
  }
  Httilde[,,1] = backCast
  
  for (t in 2:(T+1)){
    Httilde[,,t] = mu*Httilde[,,t-1] +(1-mu)*data[,,t-1]
  }
  Ht = Ht + w[k]*Httilde
}
return(Ht)
}

