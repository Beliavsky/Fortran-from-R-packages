#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]] 
#define ARMA_DONT_PRINT_ERRORS  
 
using namespace Rcpp;           
 
#ifdef _OPENMP
#include <omp.h>
//  #else
//  int mnote=1; 
#endif
 
// #define ARMA_NO_DEBUG                            
                             
//  // [[Rcpp::plugins(openmp)]]                       

// [[Rcpp::plugins("cpp11")]]  
    

#include <chrono>        
#include <iostream> 
#include <random> 
#include <cmath> 
#include <stdio.h>




// calculates the log of the rising factorial up to p
double lrfact(double n, int p)
{  double calcu=std::log(n);   // this the result for p=1
  if (p>1){
    
    for (int ii=1; ii<=(p-1); ii++)
    {calcu=calcu+std::log(n+ii);}
  }
  else if (p==0){
    calcu=std::log(1);
  }
  return calcu; 
}


void CalcuLogfac(int niter, int NIT, double n, 
                 arma::mat &alogfac, arma::mat &alogfac2, 
                 arma::mat &alfac, int nproc)
{

  
  int donde=(niter>NIT)*niter+(NIT>=niter)*NIT; 
 
  
#ifdef _OPENMP
  omp_set_num_threads(nproc);
#pragma omp parallel for   
#endif
  for (int h=0; h<=donde; h++)
  {for (int hold=0; hold<=NIT; hold++)
  { 
    alogfac(hold,h)=lrfact((n+1)*0.5+hold,h);
  }
  alogfac2(h,0)=lrfact(n*0.5,h);
    alfac(h,0)=lrfact(1,h);
  }
 }


//  Hypergeometric function that uses the rising factorials as inputs   
double ourgeoef(int h, arma::mat alogfac, arma::mat alogfac2, arma::mat alfac, double zstar, int niter=500)
{ double aux1, aux2, aux3, aux4, aux5, termo;
aux4=std::log(zstar);
double sum=1 ;
for (int s=1; s< niter; s++)
{
  aux1=alogfac(h,s);
  aux2=alogfac(0,s);
  aux3=alogfac2(s,0);
  aux5=alfac(s,0);
  termo=aux1+aux2-aux3+aux4-aux5;
  sum=sum+std::exp(termo);
  aux4=aux4+std::log(zstar);
}
return(sum);
}




//' Computes the 2F1 Hypergeometric Function
//' 
//' @name ourgeo
//'
//' @description Computes the 2F1 Hypergeometric Function
//' @usage ourgeo(a1,a2,b1,zstar,niter=500)
//' @param a1 Parameter (Real)
//' @param a2 Parameter (Real)
//' @param b1 Parameter (Real)
//' @param zstar Primary real argument
//' @param niter The degree of approximation to truncate the hypergeometric sum. The default value is set at 500
//' @return returns the value of the hypergeometric function
//' @examples
//' ##usage of ourgeo to evaluate a 2F1 hypergeometric function
//' ourgeo(1.5,1.9,1.2,0.7)
//' @export
// [[Rcpp::export]]
double ourgeo(double a1, double a2,double b1,double zstar, int niter=500)
{ int s1, s2; double aux1, aux2, aux3, aux4, aux5, termo;
s1=-1+2*(a1>0); 
aux1=std::log(std::abs(a1));
aux2=std::log(std::abs(a2));
s2=-1+2*(a2>0);
aux3=std::log(b1);
aux4=std::log(zstar);
aux5=std::log(1);
double sum=1 ;

for (int s=1; s< niter; s++)
{
   termo=aux1+aux2-aux3+aux4-aux5;
  sum=sum+s1*s2*std::exp(termo);
  aux1=aux1+std::log(std::abs(a1+s));
  s1=s1*(-1+((a1+s)>0)*2);
  aux2=aux2+std::log(std::abs(a2+s));
  s2=s2*(-1+((a2+s)>0)*2);
  aux3=aux3+std::log((b1+s));
  aux4=aux4+std::log(zstar);
  aux5=aux5+std::log(s+1);
}

return(sum);
}



//' Compute the log likelihood for an inverse gamma stochastic volatility model
//' @importFrom Rcpp evalCpp
//' 
//' @name lik_clo
//' 
//' 
//' @description Computes the log likelihood for an inverse gamma stochastic volatility model using a closed form expression of the likelihood. The details of the computation of this closed form expression are given in Leon-Gonzalez, R., & Majoni, B. (2023). Exact Likelihood for Inverse Gamma Stochastic Volatility Models (No. 23-11). Computations in 'MAC OS' are single-threaded if 'OpenMP' is not installed.
//' @details The closed form expression is obtained for the log likelihood of a stationary inverse gamma stochastic volatility model by marginalising out the volatilities. This allows the user to obtain the maximum likelihood estimator for this non linear non Gaussian state space model. When combined with `DrawK0`, the function can in addition obtain the estimates of the smoothed volatilities using the exact smoothing distributions.
//' @usage lik_clo( Res,  b2,  n, rho,  NIT=200,  niter=200,  nproc=2,  nproc2=2)
//' @param Res Matrix of OLS residuals. Usually resulting from a call to priorvar.
//' @param b2 Level of volatility.
//' @param n Degrees of freedom.
//' @param rho The parameter for the persistence of volatility.
//' @param NIT The degree of approximation to truncate the log likelihood sum. The default value is set at 200.
//' @param niter The degree of approximation to truncate the hypergeometric sum. The default value is set at 200.
//' @param nproc The number of processors allocated to evaluating the hypergeometric function. The default value is set at 2.
//' @param nproc2 The number of processors allocated to computing the log likelihood. The default value is set at 2.
//' @return A list of 7 items. List item number 1, is the sum of the log likelihood, while the rest are constants that are useful to obtain the smoothed estimates of the volatility.
//' @examples
//' ##simulate data
//' n=150
//' dat<-data.frame(Ydep=runif(n,0.3,1.4))
//' Ydep <- as.matrix(dat, -1,ncol=ncol(dat))
//' littlerho=0.95
//' r0=1
//' rho=diag(r0)*littlerho
//' p=4
//' n=4.1
//' T=nrow(Ydep)
//' Xdep <- Ydep[p:(T-1),]
//' if (p>1){
//' for(lagi in 2:p){
//'   Xdep <- cbind(Xdep, Ydep[(p-lagi+1):(T-lagi),])
//' }
//'}
//' T=nrow(Ydep)
//'  Ydep <- as.matrix(Ydep[(p+1):T,])
//'  T=nrow(Ydep)
//' unos <- rep(1,T)
//' Xdep <- cbind(unos, Xdep)
//'##obtain residuals
//' bOLS <- solve(t(Xdep) %*% Xdep) %*% t(Xdep) %*% Ydep
//' Res= Ydep- Xdep %*% bOLS
//'  Res=Res[1:T,1]
//'  b2=solve(t(Res) %*% Res/T)*(1-rho %*% rho)/(n-2)
//'  Res=as.matrix(Res,ncol=1)
//'  
//' ##obtain log likelihood
//' LL1=lik_clo(Res,b2,n,rho)
//' LL1[1]
//' @export
// [[Rcpp::export]]
Rcpp::List lik_clo(arma::mat Res, double b2,
                     double n, double rho, int NIT=200, int niter=200, int nproc=2, int nproc2=2) 
  //  NIT is the degree of approximation
{ 
  // if (mnote==1){Rcout << "omp support not found: " << " multi-thread computing disabled "}
  int T=Res.n_rows;
  arma::mat logLik=arma::zeros(T,1); 
  arma::vec oldctil=arma::zeros(NIT+1,1);
  arma::vec newctil=arma::zeros(NIT+1,1);
  arma::mat alln=arma::zeros(NIT+1,1);
  arma::mat allctil=arma::zeros(T,NIT+1); 
 
  arma::mat allc=arma::zeros(T,NIT+1); 
  int accel=0;  // whether to accelerate or not   
  int donde=(niter>NIT)*niter+(NIT>=niter)*NIT; 
  arma::mat alogfac=arma::zeros(NIT+1,donde+1);
  arma::mat alogfac2=arma::zeros(donde+1,1);
  arma::mat alfac=arma::zeros(donde+1,1); 
  
  
#ifdef _OPENMP
  omp_set_num_threads(nproc2);
#pragma omp parallel for
#else
  // single threaded version
  nproc2=1;
#endif   
  for (int h=0; h<=donde; h++)
  {for (int hold=0; hold<=NIT; hold++)
  { 
    alogfac(hold,h)=lrfact((n+1)*0.5+hold,h);
  }
  alogfac2(h,0)=lrfact(n*0.5,h);
    alfac(h,0)=lrfact(1,h);
  }

  
  double Vinv=1-rho*rho; 
  double St, Vinvtil, deltaht, normsum, liksum;
  
  // for t=0
  double et=Res(0);
  double useme=b2*et*et;
  double l0 =-0.5*(n+1)*std::log(0.5*(Vinv+useme));
  l0 =l0 + 0.5*n*std::log(0.5*Vinv) + (std::lgamma(0.5*(n+1))-std::lgamma(0.5*n))-0.5*std::log(2*arma::datum::pi)+0.5*log(b2);
  logLik(0)=l0;
  
  // for t=1
  Vinvtil=Vinv+b2*et*et;
  double bc =(0.5*rho*rho)/(rho*rho+Vinvtil);
  Vinv=1-rho*rho/(Vinvtil+rho*rho);
  
  et=Res(1);
  double deltah2=rho*rho*1.0/(rho*rho+Vinvtil);
  double  z=deltah2/(b2*et*et+1);
  // double  wz=z/(z-1);
  double  ccc=std::lgamma(0.5*(n+1));
  double c2= (std::pow(1-z,-0.5*(n+2)))*ourgeo(-0.5,-0.5,0.5*n,z,niter);
  double c2h=std::log(c2);
  l0 = -0.5*std::log(2*arma::datum::pi)+0.5*std::log(b2)+ 0.5*(n+1)*std::log(2);
  l0 = l0 + ccc - 0.5*(n+1)*std::log(b2*et*et+1);
  normsum = (std::pow(2,0.5*n))*(std::tgamma(0.5*n))*(std::pow(1-deltah2,-0.5*(n+1)));
  logLik(1)=l0+c2h-log(normsum);
  
  // for t=2
  double S3=1/(b2*et*et+1+rho*rho);
  double deltah3=(1/(1-S3*rho*rho))*S3*((rho*rho)/(Vinvtil+rho*rho));
  double c3=(std::pow(1-deltah3,-0.5*(n+2)))*ourgeo(-0.5,-0.5,0.5*n,deltah3,niter)*std::tgamma(0.5*(n+1))*std::pow(1-S3*rho*rho,-0.5*(n+1))*std::pow(2*S3,0.5*(n+1));
  
  et=Res(2);
  Vinvtil=1+b2*et*et;
  double z3=(S3*rho*rho)/Vinvtil;
  l0 = -0.5*log(2*arma::datum::pi)+0.5*log(b2)-log(c3);
  liksum=0;
  double useful=1/std::pow(Vinvtil,0.5*(n+1))*(std::tgamma(0.5*(n+1))/std::tgamma(0.5*(n)))*(std::pow(2,0.5*(n+1))/std::pow(2,0.5*n)); 
  useful=std::log(useful);
  for (int h=0; h<=NIT; h++){
    double a,b,c, auxx, chat; 
    if (accel)
    { a=-0.5*(1+2*h); b=-0.5; c=0.5*n;
    auxx=ourgeo(a,b,c,z3,niter);
    chat=-0.5*(n+2+2*h)*std::log(1-z3)+std::log(std::abs(auxx));
    }
    else {
      auxx=ourgeoef(h, alogfac,alogfac2, alfac, z3,niter);
      chat=std::log(std::abs(auxx));
    }
    double  c2til=alogfac(0,h)-alogfac2(h,0)+h*std::log(bc)-alfac(h,0); 
    int sign0=-1+2*(auxx>0);
    double aux1= c2til+ccc+alogfac(0,h)+0.5*(n+1+2*h)*std::log(2*S3)+chat+useful;
    liksum=liksum + sign0*std::exp(aux1);
    oldctil(h)=c2til;  // This is the log
  }
  logLik(2)=l0+std::log(liksum);
  allctil.row(1)=arma::trans(oldctil);
  
  //for all t 
  
  arma::mat AllVinv=arma::zeros(T,1);   arma::mat Allzt=arma::zeros(T,1);
  arma::mat AllSt=arma::zeros(T+1,1); arma::mat Alldelta=arma::zeros(T,1);  
  St=S3;
  double zt;
  AllSt(2)=St; 
  for (int tt=3; tt<=T-1; tt++){
         St=1/(b2*et*et+1+rho*rho);
    deltaht=(1/(1-St*rho*rho))*St*((rho*rho)/(Vinvtil+rho*rho));
    et=Res(tt);
    Vinvtil=1+b2*et*et;
    zt=(St*rho*rho)/Vinvtil;
    AllSt(tt)=St; Alldelta(tt)=deltaht; AllVinv(tt)=Vinvtil; Allzt(tt)=zt; 
  }
  
  arma::mat AllGeo=arma::zeros(T,(NIT+1));
  
#ifdef _OPENMP
#pragma omp parallel for   
#endif
  for (int tt=3; tt<=T-1; tt++){
    double zt=Allzt(tt);
    for  (int h=0; h<=NIT; h++){
      AllGeo(tt,h)=std::log(ourgeoef(h, alogfac,alogfac2, alfac, zt,niter));
    }
  }
  
  
  double pSt;   
  l0 = -0.5*std::log(2*arma::datum::pi)+0.5*std::log(b2);
  for (int tt=3; tt<=T; tt++){
    double ct;
    pSt=AllSt(tt-1);
    if (tt<T){
      St=AllSt(tt); deltaht=Alldelta(tt); Vinvtil=AllVinv(tt); zt=Allzt(tt);   
    }
    liksum=0;
    ct=0;
    arma::mat allik=arma::zeros(NIT+1,1);
    arma::mat alct=arma::zeros(NIT+1,1); 
    double NITper=floor((NIT+1)/nproc2);  // iterations per processor
    int remain=NIT+1-NITper*nproc2; 
    arma::vec nitvec=arma::ones<arma::vec>(nproc2)*NITper;
    if (remain>0){nitvec.row(nproc2-1)+=remain;}  
    arma::vec limits=arma::cumsum(nitvec);  arma::vec trick=arma::zeros(1,1);  limits=arma::join_cols(trick, limits);
#ifdef _OPENMP
    omp_set_num_threads(nproc2);
#pragma omp parallel for
#else
    // single threaded version
    nproc2=1;
    limits=arma::zeros(2,1); limits(0,0)=0; limits(1,0)=NIT+1;
#endif
    for (int ii=0; ii<=nproc2-1; ii++){
      for (int h=limits(ii); h<=(limits(ii+1)-1); h++)
        {
        double auxt, chat, ctil, auxi1, auxi2, auxi3; 
        
        double scale0=alogfac(NIT,h); 
        for (int hold=0; hold<=NIT; hold++)
        {
         ctil=oldctil(hold)+ccc+alogfac(0,hold)+alogfac(hold,h)+0.5*(n+1+2*hold)*log(2*pSt)-scale0; 
          newctil(h)=newctil(h)+std::exp(ctil);
          
        }
        newctil(h)=std::log(newctil(h))-alogfac2(h,0)+h*std::log(0.5*rho*rho*pSt)-alfac(h)+scale0;
        
        if (tt<T){
          auxi1=newctil(h);
         
          auxi2=-0.5*(n+1+2*h)*std::log(1-rho*rho*St);
         
          auxi3=ccc+alogfac(0,h)+0.5*(n+1+2*h)*std::log(2*St);  
          
          alct(h)=auxi1+auxi2+auxi3;
          
          if (accel){
            chat=-0.5*(n+2+2*h)*std::log(1-zt)+std::log(ourgeo(-0.5*(1+2*h),-0.5,0.5*n,zt,niter));}
          else {
           
            chat=AllGeo(tt,h);
           
          }
          
          auxt=newctil(h)+ccc+alogfac(0,h)+0.5*(n+1+2*h)*std::log(2*St)+chat;
          allik(h)=auxt; 
          
        }    
      }
    }
    
    if (tt<T){
      double cc0=alct.max()+7;    double ll0=allik.max()+7;
      ct=arma::accu(arma::exp(alct-cc0)); liksum=arma::accu(arma::exp(allik-ll0)); 
      
      liksum=liksum/(std::pow(Vinvtil,0.5*(n+1)))*((std::tgamma(0.5*(n+1))/std::tgamma(0.5*(n)))*((std::pow(2,0.5*(n+1)))/std::pow(2,0.5*n)));
      logLik(tt)=l0+std::log(liksum)-std::log(ct)+ll0-cc0;
    }
        newctil=newctil-(newctil(0));
    allctil.row(tt-1)=arma::trans(newctil); 
    oldctil=newctil ;
    newctil=arma::zeros(NIT+1,1);
    
  }
  
  // For periods 2 and (T+1)
  et=Res(0);
  AllSt(1)=1/(b2*et*et+1);
  et=Res(T-1);
  AllSt(T)=1/(b2*et*et+1);
  
  double finalLK=arma::accu(logLik);

  return Rcpp::List::create(finalLK, logLik, AllSt, allctil, alogfac, 
                            alogfac2, alfac);
} 


//' Obtains a random draw from the exact posterior of the inverse volatilities.
//'  
//' @name DrawK0
//' @description Obtains a draw from the posterior distribution of the inverse volatilities.
//' @usage DrawK0(AllSt, allctil, alogfac, alogfac2, alfac, n, rho, b2, nproc2=2)
//' @param AllSt Some constants obtained from the evaluation of the log likelihood using the function lik_clo 
//' @param allctil Some constants obtained from the evaluation of the log likelihood using the function lik_clo 
//' @param alogfac Some constants obtained from the evaluation of the log likelihood using the function lik_clo 
//' @param alogfac2 Some constants obtained from the evaluation of the log likelihood using the function lik_clo 
//' @param alfac Some constants obtained from the evaluation of the log likelihood using the function lik_clo 
//' @param n Degrees of freedom.
//' @param rho The parameter for the persistence of volatility.
//' @param b2 Level of volatility.
//' @param nproc2 The number of processors allocated to the calculations. The default value is set at 2.
//' @return A vector with a random draw from the posterior of the inverse volatilities.
//' @examples
//' ##example using US inflation Data
//' data("US_Inf_Data")
//'  Ydep <- as.matrix(US_Inf_Data)
//'  littlerho=0.95
//' r0=1
//' rho=diag(r0)*littlerho
//' p=4
//' n=4.1
//' T=nrow(Ydep)
//' Xdep <- Ydep[p:(T-1),]
//' if (p>1){
//' for(lagi in 2:p){
//'  Xdep <- cbind(Xdep, Ydep[(p-lagi+1):(T-lagi),])
//'  }
//' }
//' T=nrow(Ydep)
//' Ydep <- as.matrix(Ydep[(p+1):T,])
//'  T=nrow(Ydep)
//' unos <- rep(1,T)
//' Xdep <- cbind(unos, Xdep)
//'
//' ##obtain residuals
//'  bOLS <- solve(t(Xdep) %*% Xdep) %*% t(Xdep) %*% Ydep
//'  Res= Ydep- Xdep %*% bOLS
//'  Res=Res[1:T,1]
//'  b2=solve(t(Res) %*% Res/T)*(1-rho %*% rho)/(n-2)
//'  Res=as.matrix(Res,ncol=1)
//'   
//' ##obtain log likelihood
//'    LL1=lik_clo(Res,b2,n,rho)
//'    
//' ##obtain smoothed estimates of volatility. First, save the constants from LL1
//'  deg=200
//'  niter=200
//'  AllSt=matrix(unlist(LL1[3]), ncol=1)
//'  allctil=matrix(unlist(LL1[4]),nrow=T, ncol=(deg+1))
//'  donde=(niter>deg)*niter+(deg>=niter)*deg 
//'  alogfac=matrix(unlist(LL1[5]),nrow=(deg+1),ncol=(donde+1))
//'  alogfac2=matrix(unlist(LL1[6]), ncol=1)
//'  alfac=matrix(unlist(LL1[7]), ncol=1)
//'         
//'  milaK=0
//'  repli=5
//'   keep0=matrix(0,nrow=repli, ncol=1)
//'         for (jj in 1:repli)
//'        {
//'           laK=DrawK0(AllSt,allctil,alogfac, alogfac2, alfac, n, rho, b2, nproc2=2)
//'           
//'            milaK=milaK+1/laK*(1/repli)
//'           keep0[jj]=mean(1/laK)/b2
//'         }
//'         ccc=1/b2
//'         fefo=(milaK[1:T])*ccc
//'         
//' ##moving average of squared residuals
//'         mRes=matrix(0,nrow=T,ncol=1)
//'            Res2=Res*Res
//'          bandi=5
//'        for (iter in 1:T)
//'         {  low=(iter-bandi)*(iter>bandi)+1*(iter<=bandi)
//'           up=(iter+bandi)*(iter<=(T-bandi))+T*(iter>(T-bandi))
//'           mRes[iter]=mean(Res2[low:up])
//'         }
//'         
//' ##plot the results
//'        plot(fefo,type="l", col = "red", xlab="Time",ylab="Volatility Means")
//'           lines(mRes, type="l", col = "blue")
//'           legend("topright", legend = c("Stochastic Volatility", "Squared Residuals"),
//'                  col = c("red", "blue"), lty = 1, cex = 0.8)
//' @export
// [[Rcpp::export]]
arma::vec DrawK0(arma::mat AllSt, arma::mat allctil, arma::mat alogfac, 
                 arma::mat alogfac2, arma::mat alfac, 
                 double n, double rho, double b2, int nproc2=2)
{ 

  int NIT=allctil.n_cols-1;
  int T = allctil.n_rows;
  arma::mat allc=arma::zeros(T, NIT+1);
  arma::mat AllW=arma::zeros(T,NIT+1); 
  arma::mat AllK=arma::zeros(T,1); 
  arma::mat  rowW=arma::zeros(1,NIT+1);
  arma::mat sand;  
  // For t=T 
  int cualt=T; double sss; 
  for (int h=0; h<=(NIT); h++){
    allc(cualt-1,h)=allctil(cualt-1,h); 
    sss=allc(cualt-1,h)+alogfac(0,h)+0.5*(n+1+2*h)*std::log(2*AllSt(cualt));
    rowW(0,h)=sss; 
  }    
  double sss0=rowW.max()+7;
  rowW=arma::exp(rowW-sss0); 
  sss=arma::accu(rowW); 
  rowW=rowW*(1.0/sss);  
  AllW.row(cualt-1)=rowW; 
  
  std::mt19937 generator((std::chrono::steady_clock::now().time_since_epoch().count()));
  std::uniform_real_distribution<double> distu(0.0,1.0);
  double drawu, Khere, fredom; int where; 
  drawu=distu(generator); 
  sand=arma::cumsum(arma::trans(rowW));  
  where=std::round(arma::accu(sand<drawu));  // This gives indices from 0 to NIT+1, which is convenient in c++
  fredom=(n+1)*0.5; 
  std::gamma_distribution<double> distgam((fredom+where),2*AllSt(cualt));
  Khere=distgam(generator);
  AllK(cualt-1)=Khere;
  
  cualt=T-1; // for t=T-1, ..., 1
#ifdef _OPENMP
  omp_set_num_threads(nproc2);
#endif
  for (int cualt=(T-1); cualt>=1; cualt--)
  {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int h=0; h<=(NIT); h++){
      double auxil1=0; double auxil2, sss1; 
      if (cualt>1){
        for (int ht=0; ht<=h; ht++)
        {
          auxil2 = allctil(cualt-1,h-ht)+ht*std::log(Khere)+ht*std::log(0.25*rho*rho)-alogfac2(ht,0)-alfac(ht,0);
          auxil1+=std::exp(auxil2);    
        }
      }
      else if (cualt==1)
      {
        auxil1=h*std::log(Khere)+h*std::log(0.25*rho*rho)-alogfac2(h,0)-alfac(h,0);
        auxil1=std::exp(auxil1);
      }
      allc(cualt-1,h)=std::log(auxil1); 
      sss1=allc(cualt-1,h)+alogfac(0,h)+0.5*(n+1+2*h)*std::log(2*AllSt(cualt));
      rowW(0,h)=sss1; 
    }
    sss0=rowW.max()+7;
    rowW=arma::exp(rowW-sss0); 
    sss=arma::accu(rowW); 
    rowW=rowW*(1.0/sss);  
    AllW.row(cualt-1)=rowW; 
    
    drawu=distu(generator); 
    sand=arma::cumsum(arma::trans(rowW));  
    where=std::round(arma::accu(sand<drawu));  // This gives indices from 0 to NpartHere-1, which is convenient in c++
    std::gamma_distribution<double> distgam((fredom+where),2*AllSt(cualt));
    Khere=distgam(generator);
    AllK(cualt-1)=Khere;
  }
  return AllK;  
}



/*** R


*/
