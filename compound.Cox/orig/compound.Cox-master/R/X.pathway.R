X.pathway=function(n,p,q1,q2,rho1=0.5,rho2=0.5){
if(q1+q2>p){warning("q1+q2>p")}
#q1=sum(beta>0) ###  the number of beta>0 ###
#q2=sum(beta<0) ###  the number of beta<0 ###
q0=p-q1-q2 ### the number of beta=0 ###

r1=1/( 1+sqrt((1-rho1)/rho1) )
r2=1/( 1+sqrt((1-rho2)/rho2) )
SD0=sqrt(3/4) ## 0.866
SD1=sqrt( (3/4)*(2*r1^2-2*r1+1) )  ## 0.612 if rho1=0.5
SD2=sqrt( (3/4)*(2*r2^2-2*r2+1) )  ## 0.612 if rho2=0.5

X=matrix(0,n,p)
for(i in 1:n){
  if(q0>0){
    X[i,(q1+q2+1):p]=runif(q0,min=-1.5,max=1.5)/SD0
  }
  if(q1>0){
    A1=runif(1,min=-1.5*r1,max=1.5*r1)
    X[i,1:q1]=( runif(q1,min=-1.5*(1-r1),max=1.5*(1-r1))+A1 )/SD1
  }
  if(q2>0){
    A2=runif(1,min=-1.5*r2,max=1.5*r2)
    X[i,(q1+1):(q1+q2)]=( runif(q2,min=-1.5*(1-r2),max=1.5*(1-r2))+A2 )/SD2
  }
}
X
}
