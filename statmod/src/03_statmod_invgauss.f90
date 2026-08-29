! Inverse-Gaussian DPQR routines derived from statmod R/invgauss.R.
module statmod_invgauss
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_finite
use r_compat, only: dp, qnorm, qgamma, qchisq, pgamma, rnorm1, runif1
implicit none
private
public :: dinvgauss, pinvgauss, qinvgauss, rinvgauss
real(dp), parameter :: sqrt2 = sqrt(2.0_dp), pi=acos(-1.0_dp)
contains

pure elemental function pnorm_tail(z,lower_tail) result(p)
real(dp),intent(in)::z
logical,intent(in)::lower_tail
real(dp)::p
if(lower_tail) then
   p=0.5_dp*erfc(-z/sqrt2)
else
   p=0.5_dp*erfc(z/sqrt2)
end if
end function pnorm_tail

pure elemental function log_pnorm_tail(z,lower_tail) result(lp)
real(dp),intent(in)::z
logical,intent(in)::lower_tail
real(dp)::lp,p,zz
p=pnorm_tail(z,lower_tail)
if(p>0.0_dp) then
   lp=log(p)
else
   if(lower_tail) then
   zz=-z
   else
   zz=z
   end if
   lp=-0.5_dp*zz*zz-log(zz)-0.5_dp*log(2.0_dp*pi)
end if
end function log_pnorm_tail

pure elemental function dinvgauss(x,mean,dispersion,log_density) result(out)
real(dp),intent(in)::x
real(dp),intent(in),optional::mean,dispersion
logical,intent(in),optional::log_density
real(dp)::out,mu,phi,xx,ld
logical::lg
mu=1.0_dp
if(present(mean)) mu=mean
phi=1.0_dp
if(present(dispersion)) phi=dispersion
lg=.false.
if(present(log_density)) lg=log_density
if(x<0.0_dp .or. mu<0.0_dp .or. phi<0.0_dp) then
   out=ieee_value(out,ieee_quiet_nan)
   return
end if
if((.not.ieee_is_finite(x) .and. x>0.0_dp) .or. x<0.0_dp .or. &
   (x==0.0_dp .and. mu>0.0_dp .and. ieee_is_finite(phi))) then
   ld=-huge(1.0_dp)
else if((x>mu .and. (mu==0.0_dp .or. phi==0.0_dp)) .or. (x>0.0_dp .and. .not.ieee_is_finite(phi))) then
   ld=-huge(1.0_dp)
else if((x==mu .and. (mu==0.0_dp .or. phi==0.0_dp)) .or. (x==0.0_dp .and. .not.ieee_is_finite(phi))) then
   ld=huge(1.0_dp)
else if(.not.ieee_is_finite(mu) .and. mu>0.0_dp) then
   ld=(-log(phi)-log(2.0_dp*pi)-3.0_dp*log(x)-1.0_dp/(phi*x))/2.0_dp
else
   xx=x/mu
   ld=(-log(phi*mu)-log(2.0_dp*pi)-3.0_dp*log(xx)-(xx-1.0_dp)**2/(phi*mu*xx))/2.0_dp-log(mu)
end if
if(lg) then
out=ld
else
out=exp(ld)
end if
end function dinvgauss

function pinvgauss(q,mean,dispersion,lower_tail,log_p) result(out)
real(dp),intent(in)::q
real(dp),intent(in),optional::mean,dispersion
logical,intent(in),optional::lower_tail,log_p
real(dp)::out,mu,phi,qq,pq,a,b,lp,cv2,z
logical::lower,lg
mu=1.0_dp
if(present(mean)) mu=mean
phi=1.0_dp
if(present(dispersion)) phi=dispersion
lower=.true.
if(present(lower_tail)) lower=lower_tail
lg=.false.
if(present(log_p)) lg=log_p
if(q<0.0_dp) then
   if(lower) then
   lp=-huge(1.0_dp)
   else
   lp=0.0_dp
   end if
else if(q==0.0_dp .and. mu>0.0_dp .and. ieee_is_finite(phi)) then
   if(lower) then
   lp=-huge(1.0_dp)
   else
   lp=0.0_dp
   end if
else if(.not.ieee_is_finite(q) .and. q>0.0_dp) then
   if(lower) then
   lp=0.0_dp
   else
   lp=-huge(1.0_dp)
   end if
else if(mu<=0.0_dp .or. phi<=0.0_dp) then
   if(q<mu) then
      if(lower) then
      lp=-huge(1.0_dp)
      else
      lp=0.0_dp
      end if
   else
      if(lower) then
      lp=0.0_dp
      else
      lp=-huge(1.0_dp)
      end if
   end if
else if(.not.ieee_is_finite(mu) .and. mu>0.0_dp) then
   z=1.0_dp/(q*phi)
   if(lower) then
      ! P(X<=q)=P(chisq1 >= 1/(q phi))
      lp=log(max(0.5_dp*erfc(sqrt(z/2.0_dp)),tiny(1.0_dp)))
   else
      lp=log(max(1.0_dp-0.5_dp*erfc(sqrt(z/2.0_dp)),tiny(1.0_dp)))
   end if
else
   cv2=mu*phi
   if(cv2<1.0e-14_dp) then
      ! Gamma moment match used by upstream for tiny CV.
      qq=q/(cv2*mu)
      a=pgamma(qq,1.0_dp/cv2,rate=1.0_dp)
      if(lower) then
      lp=log(max(a,tiny(1.0_dp)))
      else
      lp=log(max(1.0_dp-a,tiny(1.0_dp)))
      end if
   else
      qq=q/mu
      phi=cv2
      pq=sqrt(phi*qq)
      a=log_pnorm_tail((qq-1.0_dp)/pq,lower)
      b=2.0_dp/phi+log_pnorm_tail(-(qq+1.0_dp)/pq,.true.)
      if(lower) then
         lp=a+log(1.0_dp+exp(b-a))
      else
         if(b>=a) then
            lp=-huge(1.0_dp)
         else
            lp=a+log(1.0_dp-exp(b-a))
         end if
         if(qq>1.0e6_dp .and. qq/(2.0_dp*phi)>5.0e5_dp) &
            lp=1.0_dp/phi-0.5_dp*log(pi)-log(2.0_dp*phi)-1.5_dp*log(1.0_dp+qq/(2.0_dp*phi))-qq/(2.0_dp*phi)
      end if
   end if
end if
if(lg) then
out=lp
else
out=exp(lp)
end if
end function pinvgauss



function rinvgauss(mean,dispersion) result(x)
real(dp),intent(in),optional::mean,dispersion
real(dp)::x,mu,phi,y,yphi,x1
mu=1.0_dp
if(present(mean)) mu=mean
phi=1.0_dp
if(present(dispersion)) phi=dispersion
if(mu<=0.0_dp .or. phi<=0.0_dp) then
   x=ieee_value(x,ieee_quiet_nan)
   return
end if
if(.not.ieee_is_finite(mu)) then
   y=rnorm1()
   x=1.0_dp/(phi*y*y)
   return
end if
if(.not.ieee_is_finite(phi)) then
x=0.0_dp
return
end if
y=rnorm1()**2
yphi=y*phi*mu
if(yphi>5.0e5_dp) then
   x1=1.0_dp/yphi
else
   x1=1.0_dp+yphi/2.0_dp*(1.0_dp-sqrt(1.0_dp+4.0_dp/yphi))
end if
if(runif1()<1.0_dp/(1.0_dp+x1)) then
x=mu*x1
else
x=mu/x1
end if
end function rinvgauss

function qinvgauss(p,mean,dispersion,lower_tail,log_p,maxit,tol) result(q)
real(dp),intent(in)::p
real(dp),intent(in),optional::mean,dispersion,tol
logical,intent(in),optional::lower_tail,log_p
integer,intent(in),optional::maxit
real(dp)::q,mu,phi,target,lp,cv2,kappa,x,dx,dens,cdf,qt,alpha,tolerance
logical::lower,lg
integer::iter,mx
mu=1.0_dp
if(present(mean)) mu=mean
phi=1.0_dp
if(present(dispersion)) phi=dispersion
lower=.true.
if(present(lower_tail)) lower=lower_tail
lg=.false.
if(present(log_p)) lg=log_p
mx=200
if(present(maxit)) mx=maxit
tolerance=1.0e-14_dp
if(present(tol)) tolerance=tol
lp=merge(p,log(max(p,tiny(1.0_dp))),lg)
target=exp(lp)
if(mu<=0.0_dp .or. phi<=0.0_dp .or. target<0.0_dp .or. target>1.0_dp) then
   q=ieee_value(q,ieee_quiet_nan)
   return
end if
if(target==0.0_dp) then
q=merge(0.0_dp,ieee_value(q,ieee_positive_inf),lower)
return
end if
if(target==1.0_dp) then
q=merge(ieee_value(q,ieee_positive_inf),0.0_dp,lower)
return
end if
if(.not.ieee_is_finite(mu)) then
   qt=merge(1.0_dp-target,target,lower)
   q=1.0_dp/(qchisq(qt,1.0_dp)*phi)
   return
end if
cv2=phi*mu
if(cv2<1.0e-8_dp) then
   qt=merge(target,1.0_dp-target,lower)
   q=qgamma(qt,1.0_dp/cv2,rate=1.0_dp/(cv2*mu))
   return
end if
phi=cv2
kappa=1.5_dp*phi
x=sqrt(1.0_dp+kappa*kappa)-kappa
if(kappa>1.0e3_dp) then
qt=1.0_dp/(2.0_dp*kappa)
x=qt*(1.0_dp-qt*qt)
end if
if((lower .and. lp < -11.51_dp) .or. (.not.lower .and. lp > -1.0e-5_dp)) then
   qt=qnorm(target,lower_tail=lower)
   x=1.0_dp/(phi*qt*qt)
end if
if((lower .and. lp > -1.0e-5_dp) .or. (.not.lower .and. lp < -11.51_dp)) then
   alpha=1.0_dp/phi
   qt=merge(target,1.0_dp-target,lower)
   x=max(x,qgamma(qt,alpha,rate=alpha))
end if
do iter=1,mx
   cdf=pinvgauss(x,dispersion=phi,lower_tail=lower)
   dens=dinvgauss(x,dispersion=phi)
   if(lower) then
   dx=(target-cdf)/dens
   else
   dx=(target-cdf)/(-dens)
   end if
   if(.not.ieee_is_finite(dx)) exit
   if(lower) then
   x=max(tiny(1.0_dp),x+dx)
   else
   x=max(tiny(1.0_dp),x+dx)
   end if
   if(abs(dx)/max(x,1.0_dp)<tolerance) exit
end do
q=mu*x
end function qinvgauss

end module statmod_invgauss
