module ld_distributions_extra
use ld_kinds, only: dp, pi
use ld_random, only: rand_uniform, rand_normal, rand_gamma, rand_chisq
use ld_distributions, only: normal_logpdf
implicit none
private
public :: dalaplace, palaplace, qalaplace, ralaplace
public :: dallaplace, pallaplace, qallaplace, rallaplace
public :: dllaplace, pllaplace, qllaplace, rllaplace
public :: dslaplace, pslaplace, qslaplace, rslaplace
public :: dsdlaplace, psdlaplace, qsdlaplace
public :: dpe, rpe, dhalft, rhalft, dinvchisq, rinvchisq, dinvgaussian, rinvgaussian
contains
pure function dalaplace(x,location,scale,kappa,log_density) result(v)
   real(dp),intent(in)::x,location,scale,kappa
   logical,intent(in),optional::log_density
   real(dp)::v,k,lc
   logical::lg
   lg=.false.; if(present(log_density)) lg=log_density
   if(scale<=0.0_dp .or. kappa<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   lc=0.5_dp*log(2.0_dp)-log(scale)+log(kappa)-log(1.0_dp+kappa*kappa)
   k=kappa; if(x<location) k=1.0_dp/kappa
   v=lc-(sqrt(2.0_dp)/scale)*abs(x-location)*k
   if(.not.lg) v=exp(v)
end function dalaplace
pure function palaplace(q,location,scale,kappa) result(p)
   real(dp),intent(in)::q,location,scale,kappa; real(dp)::p,k,e,t
   k=kappa; if(q<location) k=1.0_dp/kappa
   e=exp(-(sqrt(2.0_dp)/scale)*abs(q-location)*k); t=e/(1.0_dp+kappa*kappa)
   if(q<location) then; p=kappa*kappa*t; else; p=1.0_dp-t; end if
end function palaplace
pure function qalaplace(p,location,scale,kappa) result(q)
   real(dp),intent(in)::p,location,scale,kappa; real(dp)::q,t
   if(p<=0.0_dp) then; q=-huge(1.0_dp); return; end if
   if(p>=1.0_dp) then; q=huge(1.0_dp); return; end if
   t=kappa*kappa/(1.0_dp+kappa*kappa)
   if(p<=t) then
      q=location+scale*kappa*log(p/t)/sqrt(2.0_dp)
   else
      q=location-(scale/kappa)*(log(1.0_dp+kappa*kappa)+log(1.0_dp-p))/sqrt(2.0_dp)
   end if
end function qalaplace
function ralaplace(location,scale,kappa) result(x)
   real(dp),intent(in)::location,scale,kappa; real(dp)::x
   x=location+scale*log(rand_uniform()**kappa/rand_uniform()**(1.0_dp/kappa))/sqrt(2.0_dp)
end function ralaplace

pure function dallaplace(x,location,scale,kappa,log_density) result(v)
   real(dp),intent(in)::x,location,scale,kappa; logical,intent(in),optional::log_density
   real(dp)::v,a,b,d,e; logical::lg
   lg=.false.; if(present(log_density)) lg=log_density
   if(x<=0.0_dp .or. scale<=0.0_dp .or. kappa<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   a=sqrt(2.0_dp)*kappa/scale; b=sqrt(2.0_dp)/(scale*kappa); d=exp(location)
   if(x>=d) then; e=(b-1.0_dp)*(log(x)-location); else; e=-(a+1.0_dp)*(log(x)-location); end if
   v=-location+log(a)+log(b)-log(a+b)+e; if(.not.lg) v=exp(v)
end function dallaplace
pure function pallaplace(q,location,scale,kappa) result(p)
   real(dp),intent(in)::q,location,scale,kappa; real(dp)::p,a,b,d,t
   if(q<=0.0_dp) then; p=0.0_dp; return; end if
   a=sqrt(2.0_dp)*kappa/scale; b=sqrt(2.0_dp)/(scale*kappa); d=exp(location); t=a+b
   if(q<d) then; p=(a/t)*(q/d)**b; else; p=1.0_dp-(b/t)*(d/q)**a; end if
end function pallaplace
pure function qallaplace(p,location,scale,kappa) result(q)
   real(dp),intent(in)::p,location,scale,kappa; real(dp)::q,a,b,d,t
   if(p<=0.0_dp) then; q=0.0_dp; return; end if
   if(p>=1.0_dp) then; q=huge(1.0_dp); return; end if
   a=sqrt(2.0_dp)*kappa/scale; b=sqrt(2.0_dp)/(scale*kappa); d=exp(location); t=a+b
   if(p<=a/t) then; q=d*(p*t/a)**(1.0_dp/b); else; q=d*((1.0_dp-p)*t/b)**(-1.0_dp/a); end if
end function qallaplace
function rallaplace(location,scale,kappa) result(x)
   real(dp),intent(in)::location,scale,kappa; real(dp)::x
   x=exp(location)*(rand_uniform()**kappa/rand_uniform()**(1.0_dp/kappa))**(scale/sqrt(2.0_dp))
end function rallaplace

pure function dllaplace(x,location,scale,log_density) result(v)
   real(dp),intent(in)::x,location,scale; logical,intent(in),optional::log_density; real(dp)::v
   v=dallaplace(x,location,scale,1.0_dp,log_density)
end function dllaplace
pure function pllaplace(q,location,scale) result(p)
   real(dp),intent(in)::q,location,scale; real(dp)::p
   p=pallaplace(q,location,scale,1.0_dp)
end function pllaplace
pure function qllaplace(p,location,scale) result(q)
   real(dp),intent(in)::p,location,scale; real(dp)::q
   q=qallaplace(p,location,scale,1.0_dp)
end function qllaplace
function rllaplace(location,scale) result(x)
   real(dp),intent(in)::location,scale; real(dp)::x
   x=rallaplace(location,scale,1.0_dp)
end function rllaplace

pure function dslaplace(x,mu,alpha,beta,log_density) result(v)
   real(dp),intent(in)::x,mu,alpha,beta; logical,intent(in),optional::log_density; real(dp)::v; logical::lg
   lg=.false.; if(present(log_density)) lg=log_density
   if(alpha<=0.0_dp .or. beta<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   if(x<=mu) then; v=-log(alpha+beta)+(x-mu)/alpha; else; v=-log(alpha+beta)+(mu-x)/beta; end if
   if(.not.lg) v=exp(v)
end function dslaplace
pure function pslaplace(q,mu,alpha,beta) result(p)
   real(dp),intent(in)::q,mu,alpha,beta; real(dp)::p
   if(q<mu) then; p=alpha/(alpha+beta)*exp((q-mu)/alpha); else; p=1.0_dp-beta/(alpha+beta)*exp((mu-q)/beta); end if
end function pslaplace
pure function qslaplace(p,mu,alpha,beta) result(q)
   real(dp),intent(in)::p,mu,alpha,beta; real(dp)::q,t
   t=alpha/(alpha+beta)
   if(p<t) then; q=alpha*log(p*(alpha+beta)/alpha)+mu; else; q=mu-beta*log((alpha+beta)*(1.0_dp-p)/beta); end if
end function qslaplace
function rslaplace(mu,alpha,beta) result(x)
   real(dp),intent(in)::mu,alpha,beta; real(dp)::x,e
   e=-log(rand_uniform())
   if(rand_uniform()<=alpha/(alpha+beta)) then; x=mu-alpha*e; else; x=mu+beta*e; end if
end function rslaplace

pure function dsdlaplace(x,p,q,log_density) result(v)
   integer,intent(in)::x; real(dp),intent(in)::p,q; logical,intent(in),optional::log_density; real(dp)::v; logical::lg
   lg=.false.; if(present(log_density)) lg=log_density
   if(p<0.0_dp .or. p>=1.0_dp .or. q<0.0_dp .or. q>=1.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   if(x>=0) then; v=log(1.0_dp-p)+log(1.0_dp-q)-log(1.0_dp-p*q)-real(x,dp)*log(max(p,tiny(1.0_dp)))
   else; v=log(1.0_dp-p)+log(1.0_dp-q)-log(1.0_dp-p*q)-real(abs(x),dp)*log(max(q,tiny(1.0_dp))); end if
   ! Upstream code uses subtraction in log form; preserve the implemented package formula exactly below.
   if(x>=0) v=log(1.0_dp-p)+log(1.0_dp-q)-log(1.0_dp-p*q)+real(x,dp)*log(max(p,tiny(1.0_dp)))
   if(x<0) v=log(1.0_dp-p)+log(1.0_dp-q)-log(1.0_dp-p*q)+real(abs(x),dp)*log(max(q,tiny(1.0_dp)))
   if(.not.lg) v=exp(v)
end function dsdlaplace
pure function psdlaplace(x,p,q) result(pr)
   integer,intent(in)::x; real(dp),intent(in)::p,q; real(dp)::pr
   if(x>=0) then; pr=1.0_dp-(1.0_dp-q)*p**real(x+1,dp)/(1.0_dp-p*q)
   else; pr=(1.0_dp-p)*q**real(-x,dp)/(1.0_dp-p*q); end if
end function psdlaplace
function qsdlaplace(prob,p,q) result(x)
   real(dp),intent(in)::prob,p,q; integer::x
   x=0
   if(prob>=psdlaplace(0,p,q)) then
      do while(prob>=psdlaplace(x,p,q)); x=x+1; if(x>100000) exit; end do
   else
      do while(prob<psdlaplace(x,p,q)); x=x-1; if(x< -100000) exit; end do
      x=x+1
   end if
end function qsdlaplace

pure function dpe(x,mu,sigma,kappa,log_density) result(v)
   real(dp),intent(in)::x,mu,sigma,kappa; logical,intent(in),optional::log_density; real(dp)::v; logical::lg
   lg=.false.; if(present(log_density)) lg=log_density
   if(sigma<=0.0_dp .or. kappa<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   v=-log(2.0_dp*kappa**(1.0_dp/kappa)*gamma(1.0_dp+1.0_dp/kappa)*sigma)&
      -abs(x-mu)**kappa/(kappa*sigma**kappa)
   if(.not.lg) v=exp(v)
end function dpe
function rpe(mu,sigma,kappa) result(x)
   real(dp),intent(in)::mu,sigma,kappa; real(dp)::x,z
   z=rand_gamma(1.0_dp/kappa,kappa)**(1.0_dp/kappa); if(rand_uniform()<0.5_dp) z=-z; x=mu+sigma*z
end function rpe

pure function dhalft(x,scale,nu,log_density) result(v)
   real(dp),intent(in)::x,scale,nu; logical,intent(in),optional::log_density; real(dp)::v; logical::lg
   lg=.false.; if(present(log_density)) lg=log_density
   if(x<0.0_dp .or. scale<=0.0_dp .or. nu<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   v=log(2.0_dp)+log_gamma(0.5_dp*(nu+1.0_dp))-log_gamma(0.5_dp*nu)-0.5_dp*log(pi*nu)-log(scale)&
      -0.5_dp*(nu+1.0_dp)*log(1.0_dp+(x/scale)**2/nu)
   if(.not.lg) v=exp(v)
end function dhalft
function rhalft(scale,nu) result(x)
   real(dp),intent(in)::scale,nu; real(dp)::x
   x=abs(scale*rand_normal()/sqrt(rand_chisq(nu)/nu))
end function rhalft

pure function dinvchisq(x,df,scale,log_density) result(v)
   real(dp),intent(in)::x,df,scale; logical,intent(in),optional::log_density; real(dp)::v,a,b; logical::lg
   lg=.false.; if(present(log_density)) lg=log_density; a=0.5_dp*df; b=0.5_dp*df*scale
   if(x<=0.0_dp .or. df<=0.0_dp .or. scale<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   v=a*log(b)-log_gamma(a)-(a+1.0_dp)*log(x)-b/x; if(.not.lg) v=exp(v)
end function dinvchisq
function rinvchisq(df,scale) result(x)
   real(dp),intent(in)::df,scale; real(dp)::x
   x=df*scale/rand_chisq(df)
end function rinvchisq

pure function dinvgaussian(x,mu,lambda,log_density) result(v)
   real(dp),intent(in)::x,mu,lambda; logical,intent(in),optional::log_density; real(dp)::v; logical::lg
   lg=.false.; if(present(log_density)) lg=log_density
   if(x<=0.0_dp .or. mu<=0.0_dp .or. lambda<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   v=0.5_dp*(log(lambda)-log(2.0_dp*pi)-3.0_dp*log(x))-lambda*(x-mu)**2/(2.0_dp*mu*mu*x)
   if(.not.lg) v=exp(v)
end function dinvgaussian
function rinvgaussian(mu,lambda) result(x)
   real(dp),intent(in)::mu,lambda; real(dp)::x,y,z,u
   y=rand_normal()**2
   z=mu+mu*mu*y/(2.0_dp*lambda)-mu/(2.0_dp*lambda)*sqrt(4.0_dp*mu*lambda*y+mu*mu*y*y)
   u=rand_uniform(); if(u<=mu/(mu+z)) then; x=z; else; x=mu*mu/z; end if
end function rinvgaussian
end module ld_distributions_extra
