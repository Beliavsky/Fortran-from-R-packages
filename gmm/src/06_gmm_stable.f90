! SPDX-License-Identifier: GPL-2.0-or-later
module gmm_stable
use r_compat, only: dp
implicit none
private
public :: char_stable
contains
pure function char_stable(theta,tau,pm) result(phi)
real(dp),intent(in)::theta(4),tau(:)
integer,intent(in),optional::pm
complex(dp)::phi(size(tau))
real(dp)::a,b,g,d,re,im,t,ang,pi
integer::i,p
pi=acos(-1.0_dp)
a=theta(1)
b=theta(2)
g=theta(3)
d=theta(4)
p=0
if(present(pm)) p=pm
do i=1,size(tau)
   t=tau(i)
   if(g==0.0_dp) then
      phi(i)=exp(cmplx(0.0_dp,d*t,dp))
      cycle
   end if
   if(p==0) then
      if(a==1.0_dp) then
         re=-g*abs(t)
         im=d*t
         if(t/=0.0_dp) im=im+re*(2.0_dp/pi)*b*sign(1.0_dp,t)*log(g*abs(t))
      else
         ang=tan(pi*a/2.0_dp)
         re=-(g*abs(t))**a
         im=d*t
         if(t/=0.0_dp) im=im+re*b*ang*sign(1.0_dp,t)*((g*abs(t))**(1.0_dp-a)-1.0_dp)
      end if
   else
      if(a==1.0_dp) then
         re=-g*abs(t)
         im=d*t
         if(t/=0.0_dp) im=im+re*b*(2.0_dp/pi)*sign(1.0_dp,t)*log(abs(t))
      else
         ang=tan(pi*a/2.0_dp)
         re=-(g*abs(t))**a
         im=re*(-ang*b*sign(1.0_dp,t))+d*t
      end if
   end if
   phi(i)=exp(cmplx(re,im,dp))
end do
end function char_stable
end module gmm_stable
