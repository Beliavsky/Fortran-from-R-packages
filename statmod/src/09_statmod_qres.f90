module statmod_qres
use r_compat, only: dp, pbinom, ppois, pgamma, pbeta, qnorm, runif1
use statmod_invgauss, only: pinvgauss
use tweedie_distribution_mod, only: ptweedie
implicit none
private
public :: qres_binom, qres_pois, qres_gamma, qres_invgauss, qres_nbinom, qres_tweedie, qres_default
contains

subroutine qres_binom(y,p,n,resid)
real(dp),intent(in)::y(:),p(:),n(:)
real(dp),allocatable,intent(out)::resid(:)
real(dp)::a,b,u,yy
integer::i,nn
allocate(resid(size(y)))
do i=1,size(y)
   nn=nint(n(i))
   yy=n(i)*y(i)
   a=pbinom(yy-1.0_dp,nn,p(i))
   b=pbinom(yy,nn,p(i))
   u=a+(b-a)*runif1()
   resid(i)=qnorm(u)
end do
end subroutine

subroutine qres_pois(y,mu,resid)
real(dp),intent(in)::y(:),mu(:)
real(dp),allocatable,intent(out)::resid(:)
real(dp)::a,b,u
integer::i
allocate(resid(size(y)))
do i=1,size(y)
a=ppois(y(i)-1,lambda=mu(i))
b=ppois(y(i),lambda=mu(i))
u=a+(b-a)*runif1()
resid(i)=qnorm(u)
end do
end subroutine

subroutine qres_gamma(y,mu,df,weights,dispersion,resid)
real(dp),intent(in)::y(:),mu(:)
integer,intent(in)::df
real(dp),intent(in),optional::weights(:),dispersion
real(dp),allocatable,intent(out)::resid(:)
real(dp),allocatable::w(:)
real(dp)::phi,u
integer::i
allocate(w(size(y)))
w=1
if(present(weights))w=weights
if(present(dispersion))then
phi=dispersion
else
phi=sum(w*((y-mu)/mu)**2)/real(df,dp)
end if
allocate(resid(size(y)))
do i=1,size(y)
u=pgamma(w(i)*y(i)/mu(i)/phi,shape=w(i)/phi)
resid(i)=qnorm(u)
end do
end subroutine

subroutine qres_invgauss(y,mu,df,weights,dispersion,resid)
real(dp),intent(in)::y(:),mu(:)
integer,intent(in)::df
real(dp),intent(in),optional::weights(:),dispersion
real(dp),allocatable,intent(out)::resid(:)
real(dp),allocatable::w(:)
real(dp)::phi,u
integer::i
allocate(w(size(y)))
w=1
if(present(weights))w=weights
if(present(dispersion))then
phi=dispersion
else
phi=sum(w*(y-mu)**2/(mu*mu*y))/real(df,dp)
end if
allocate(resid(size(y)))
do i=1,size(y)
   if(y(i)<mu(i))then
   u=pinvgauss(y(i),mean=mu(i),dispersion=phi)
   resid(i)=qnorm(u)
   else if(y(i)>mu(i))then
   u=pinvgauss(y(i),mean=mu(i),dispersion=phi,lower_tail=.false.)
   resid(i)=qnorm(u,lower_tail=.false.)
   else
   resid(i)=0
   end if
end do
end subroutine

subroutine qres_nbinom(y,mu,sizepar,resid)
real(dp),intent(in)::y(:),mu(:),sizepar
real(dp),allocatable,intent(out)::resid(:)
real(dp)::p,a,b,u
integer::i
allocate(resid(size(y)))
do i=1,size(y)
   p=sizepar/(mu(i)+sizepar)
   if(y(i)>0)then
   a=pbeta(p,sizepar,max(y(i),1.0_dp))
   else
   a=0
   end if
   b=pbeta(p,sizepar,y(i)+1.0_dp)
   u=a+(b-a)*runif1()
   resid(i)=qnorm(u)
end do
end subroutine

subroutine qres_tweedie(y,mu,power,df,weights,dispersion,resid)
real(dp),intent(in)::y(:),mu(:),power
integer,intent(in)::df
real(dp),intent(in),optional::weights(:),dispersion
real(dp),allocatable,intent(out)::resid(:)
real(dp),allocatable::w(:)
real(dp)::phi,u
integer::i
allocate(w(size(y)))
w=1
if(present(weights))w=weights
if(present(dispersion))then
phi=dispersion
else
phi=sum(w*(y-mu)**2/mu**power)/real(df,dp)
end if
allocate(resid(size(y)))
do i=1,size(y)
   u=ptweedie(y(i),mu(i),phi/w(i),power)
   if(power>1.and.power<2.and.y(i)==0)u=runif1()*u
   resid(i)=qnorm(u)
end do
end subroutine

subroutine qres_default(dev_resid,weights,df,dispersion,resid)
real(dp),intent(in)::dev_resid(:),weights(:)
integer,intent(in)::df
real(dp),intent(in),optional::dispersion
real(dp),allocatable,intent(out)::resid(:)
real(dp)::phi
if(present(dispersion))then
phi=dispersion
else
if(df>0)then
phi=sum(weights*dev_resid**2)/df
else
phi=1
end if
end if
resid=dev_resid/sqrt(phi)
end subroutine

end module statmod_qres
