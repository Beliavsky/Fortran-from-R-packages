! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
module robustbase_psi
   use robustbase_kinds, only: dp
   implicit none
   private
   public :: huber_psi, huber_rho, huber_weight, hampel_psi, hampel_rho, hampel_weight, &
             tukey_psi, tukey_rho, tukey_weight, hard_rejection_weight, smooth_weight, &
             welsh_psi, welsh_rho, welsh_weight, welsh_psi_derivative, &
             optimal_psi, optimal_rho, optimal_weight, optimal_psi_derivative, &
             ggw_psi, ggw_rho, ggw_weight, ggw_psi_derivative, &
             lqq_psi, lqq_rho, lqq_weight, lqq_psi_derivative, lqq_psi_second_derivative
contains
   elemental function huber_psi(x,k) result(y)
      real(dp),intent(in)::x,k
      real(dp)::y
      y=max(-k,min(k,x))
   end function huber_psi

   elemental function huber_rho(x,k) result(y)
      real(dp),intent(in)::x,k
      real(dp)::y
      if(abs(x)<=k) then
         y=0.5_dp*x*x
      else
         y=k*abs(x)-0.5_dp*k*k
      end if
   end function huber_rho

   elemental function huber_weight(x,k) result(y)
      real(dp),intent(in)::x,k
      real(dp)::y
      if(abs(x)<=k .or. abs(x)<=tiny(1.0_dp)) then
         y=1.0_dp
      else
         y=k/abs(x)
      end if
   end function huber_weight

   elemental function hampel_psi(x,a,b,c) result(y)
      real(dp),intent(in)::x,a,b,c
      real(dp)::y,z
      z=abs(x)
      if(z<=a) then
         y=x
      else if(z<=b) then
         y=sign(a,x)
      else if(z<=c) then
         y=sign(a*(c-z)/(c-b),x)
      else
         y=0.0_dp
      end if
   end function hampel_psi

   elemental function hampel_rho(x,a,b,c) result(y)
      real(dp),intent(in)::x,a,b,c
      real(dp)::y,z
      z=abs(x)
      if(z<=a) then
         y=0.5_dp*z*z
      else if(z<=b) then
         y=a*z-0.5_dp*a*a
      else if(z<=c) then
         y=a*b-0.5_dp*a*a+a*(c*(z-b)-0.5_dp*(z*z-b*b))/(c-b)
      else
         y=a*b-0.5_dp*a*a+0.5_dp*a*(c-b)
      end if
   end function hampel_rho

   elemental function hampel_weight(x,a,b,c) result(y)
      real(dp),intent(in)::x,a,b,c
      real(dp)::y
      if(abs(x)<=tiny(1.0_dp)) then
         y=1.0_dp
      else
         y=hampel_psi(x,a,b,c)/x
      end if
   end function hampel_weight

   elemental function tukey_psi(x,c) result(y)
      real(dp),intent(in)::x,c
      real(dp)::y,u
      u=x/c
      if(abs(u)<1.0_dp) then
         y=x*(1.0_dp-u*u)**2
      else
         y=0.0_dp
      end if
   end function tukey_psi

   elemental function tukey_rho(x,c) result(y)
      real(dp),intent(in)::x,c
      real(dp)::y,u
      u=x/c
      if(abs(u)<1.0_dp) then
         y=c*c/6.0_dp*(1.0_dp-(1.0_dp-u*u)**3)
      else
         y=c*c/6.0_dp
      end if
   end function tukey_rho

   elemental function tukey_weight(x,c) result(y)
      real(dp),intent(in)::x,c
      real(dp)::y,u
      u=x/c
      if(abs(u)<1.0_dp) then
         y=(1.0_dp-u*u)**2
      else
         y=0.0_dp
      end if
   end function tukey_weight

   elemental function welsh_rho(x,c) result(y)
      real(dp),intent(in)::x,c
      real(dp)::y,u
      u=x/c
      y=1.0_dp-exp(-0.5_dp*u*u)
   end function welsh_rho

   elemental function welsh_psi(x,c) result(y)
      real(dp),intent(in)::x,c
      real(dp)::y,u,e
      u=x/c
      if(abs(u)>37.7_dp) then
         y=0.0_dp
      else
         e=exp(-0.5_dp*u*u)
         y=x*e
      end if
   end function welsh_psi

   elemental function welsh_psi_derivative(x,c) result(y)
      real(dp),intent(in)::x,c
      real(dp)::y,u
      u=x/c
      if(abs(u)>37.7_dp) then
         y=0.0_dp
      else
         y=exp(-0.5_dp*u*u)*(1.0_dp-u*u)
      end if
   end function welsh_psi_derivative

   elemental function welsh_weight(x,c) result(y)
      real(dp),intent(in)::x,c
      real(dp)::y,u
      u=x/c
      y=exp(-0.5_dp*u*u)
   end function welsh_weight

   elemental function optimal_rho(x,c) result(y)
      real(dp),intent(in)::x,c
      real(dp)::y,u,z
      real(dp),parameter::r1=-1.944_dp/2.0_dp,r2=1.728_dp/4.0_dp,r3=-0.312_dp/6.0_dp,r4=0.016_dp/8.0_dp
      u=x/c
      z=abs(u)
      if(z>3.0_dp) then
         y=1.0_dp
      else if(z>2.0_dp) then
         z=z*z
         y=(z*(r1+z*(r2+z*(r3+z*r4)))+1.792_dp)/3.25_dp
      else
         y=u*u/6.5_dp
      end if
   end function optimal_rho

   elemental function optimal_psi(x,c) result(y)
      real(dp),intent(in)::x,c
      real(dp)::y,u,z,u2,poly
      real(dp),parameter::r1=-1.944_dp,r2=1.728_dp,r3=-0.312_dp,r4=0.016_dp
      u=x/c
      z=abs(u)
      if(z>3.0_dp) then
         y=0.0_dp
      else if(z>2.0_dp) then
         u2=u*u
         poly=c*((((r4*u2+r3)*u2+r2)*u2+r1)*u)
         if(u>0.0_dp) then
            y=max(0.0_dp,poly)
         else
            y=-abs(poly)
         end if
      else
         y=x
      end if
   end function optimal_psi

   elemental function optimal_psi_derivative(x,c) result(y)
      real(dp),intent(in)::x,c
      real(dp)::y,u,z,z2
      real(dp),parameter::r1=-1.944_dp,r2=1.728_dp,r3=-0.312_dp,r4=0.016_dp
      u=x/c
      z=abs(u)
      if(z>3.0_dp) then
         y=0.0_dp
      else if(z>2.0_dp) then
         z2=z*z
         y=r1+z2*(3.0_dp*r2+z2*(5.0_dp*r3+z2*7.0_dp*r4))
      else
         y=1.0_dp
      end if
   end function optimal_psi_derivative

   elemental function optimal_weight(x,c) result(y)
      real(dp),intent(in)::x,c
      real(dp)::y,u,z,z2
      real(dp),parameter::r1=-1.944_dp,r2=1.728_dp,r3=-0.312_dp,r4=0.016_dp
      u=x/c
      z=abs(u)
      if(z>3.0_dp) then
         y=0.0_dp
      else if(z>2.0_dp) then
         z2=z*z
         y=max(0.0_dp,r1+z2*(r2+z2*(r3+z2*r4)))
      else
         y=1.0_dp
      end if
   end function optimal_weight

   elemental function ggw_psi(x,a,b,c) result(y)
      real(dp),intent(in)::x,a,b,c
      real(dp)::y,z,e
      z=abs(x)
      if(z<c) then
         y=x
      else
         e=-(z-c)**b/(2.0_dp*a)
         if(e<-708.4_dp) then
            y=0.0_dp
         else
            y=x*exp(e)
         end if
      end if
   end function ggw_psi

   elemental function ggw_psi_derivative(x,a,b,c) result(y)
      real(dp),intent(in)::x,a,b,c
      real(dp)::y,z,e
      z=abs(x)
      if(z<c) then
         y=1.0_dp
      else
         e=-(z-c)**b/(2.0_dp*a)
         if(e<-708.4_dp) then
            y=0.0_dp
         else if(abs(z-c)<=epsilon(1.0_dp)*max(1.0_dp,c) .and. b<1.0_dp) then
            y=1.0_dp
         else
            y=exp(e)*(1.0_dp-b*z*max(z-c,0.0_dp)**(b-1.0_dp)/(2.0_dp*a))
         end if
      end if
   end function ggw_psi_derivative

   elemental function ggw_weight(x,a,b,c) result(y)
      real(dp),intent(in)::x,a,b,c
      real(dp)::y,z
      z=abs(x)
      if(z<c) then
         y=1.0_dp
      else
         y=exp(-(z-c)**b/(2.0_dp*a))
      end if
   end function ggw_weight

   function ggw_rho(x,a,b,c,rho_inf,n_intervals) result(y)
      real(dp),intent(in)::x,a,b,c
      real(dp),intent(in),optional::rho_inf
      integer,intent(in),optional::n_intervals
      real(dp)::y,norm,h,z,sumv,upper
      integer::n,i
      n=400
      if(present(n_intervals))n=max(20,n_intervals)
      if(mod(n,2)/=0)n=n+1
      upper=abs(x)
      if(upper<=tiny(1.0_dp))then
         y=0.0_dp
         return
      end if
      h=upper/real(n,dp)
      sumv=ggw_psi(0.0_dp,a,b,c)+ggw_psi(upper,a,b,c)
      do i=1,n-1
         z=h*real(i,dp)
         sumv=sumv+merge(4.0_dp,2.0_dp,mod(i,2)==1)*ggw_psi(z,a,b,c)
      end do
      y=h*sumv/3.0_dp
      if(present(rho_inf))then
         norm=max(rho_inf,tiny(1.0_dp))
      else
         norm=ggw_rho_infinity(a,b,c)
      end if
      y=min(1.0_dp,max(0.0_dp,y/norm))
   end function ggw_rho

   function ggw_rho_infinity(a,b,c) result(value)
      real(dp),intent(in)::a,b,c
      real(dp)::value,upper,h,z,sumv
      integer,parameter::n=1600
      integer::i
      upper=c+max(20.0_dp,(2.0_dp*a*36.0_dp)**(1.0_dp/b))
      h=upper/real(n,dp)
      sumv=ggw_psi(0.0_dp,a,b,c)+ggw_psi(upper,a,b,c)
      do i=1,n-1
         z=h*real(i,dp)
         sumv=sumv+merge(4.0_dp,2.0_dp,mod(i,2)==1)*ggw_psi(z,a,b,c)
      end do
      value=max(h*sumv/3.0_dp,tiny(1.0_dp))
   end function ggw_rho_infinity

   elemental function lqq_psi_derivative(x,b,c,s) result(y)
      real(dp),intent(in)::x,b,c,s
      real(dp)::y,z,bc,a,one_minus_s
      z=abs(x)
      if(z<=c)then
         y=1.0_dp
      else
         bc=b+c
         if(z<=bc)then
            y=1.0_dp-s*(z-c)/b
         else
            one_minus_s=1.0_dp-s
            a=(b*s-2.0_dp*bc)/one_minus_s
            if(z<bc+a)then
               y=-one_minus_s*((z-bc)/a-1.0_dp)
            else
               y=0.0_dp
            end if
         end if
      end if
   end function lqq_psi_derivative

   elemental function lqq_psi_second_derivative(x,b,c,s) result(y)
      real(dp),intent(in)::x,b,c,s
      real(dp)::y,z,bc,a,one_minus_s,sx
      z=abs(x)
      sx=merge(-1.0_dp,1.0_dp,x<0.0_dp)
      if(z<=c)then
         y=0.0_dp
      else
         bc=b+c
         if(z<=bc)then
            y=-sx*s/b
         else
            one_minus_s=1.0_dp-s
            a=(b*s-2.0_dp*bc)/one_minus_s
            if(z<bc+a)then
               y=-sx*one_minus_s/a
            else
               y=0.0_dp
            end if
         end if
      end if
   end function lqq_psi_second_derivative

   elemental function lqq_psi(x,b,c,s) result(y)
      real(dp),intent(in)::x,b,c,s
      real(dp)::y,z,bc,s5,s6,sgn
      z=abs(x)
      sgn=merge(-1.0_dp,1.0_dp,x<0.0_dp)
      if(z<=c)then
         y=x
      else
         bc=b+c
         if(z<=bc)then
            y=sgn*(z-s*(z-c)**2/(2.0_dp*b))
         else
            s5=s-1.0_dp
            s6=-2.0_dp*bc+b*s
            if(z<bc-s6/s5)then
               y=sgn*(-s6/2.0_dp-s5*s5/s6*((z-bc)**2/2.0_dp+s6/s5*(z-bc)))
            else
               y=0.0_dp
            end if
         end if
      end if
   end function lqq_psi

   elemental function lqq_rho(x,b,c,s) result(y)
      real(dp),intent(in)::x,b,c,s
      real(dp)::y,z,bc,s0,s5,s6,s7,denom
      z=abs(x)
      bc=b+c
      denom=s*c*(3.0_dp*c+2.0_dp*b)+bc*bc
      if(z<=c)then
         y=(3.0_dp*s-3.0_dp)*x*x/denom
      else if(z<=bc)then
         s0=z-c
         y=(6.0_dp*s-6.0_dp)*(x*x/2.0_dp-s*s0**3/(6.0_dp*b))/denom
      else
         s5=s-1.0_dp
         s6=-2.0_dp*bc+b*s
         if(z<bc-s6/s5)then
            s7=z-bc
            y=6.0_dp*s5*(bc*bc/2.0_dp-s*b*b/6.0_dp-s7/2.0_dp*(s6+s7*(s5+s7*s5*s5/(3.0_dp*s6))))/denom
         else
            y=1.0_dp
         end if
      end if
   end function lqq_rho

   elemental function lqq_weight(x,b,c,s) result(y)
      real(dp),intent(in)::x,b,c,s
      real(dp)::y,z,bc,s0,s5,s6,s7
      z=abs(x)
      if(z<=c .or. z<=tiny(1.0_dp))then
         y=1.0_dp
      else
         bc=b+c
         if(z<=bc)then
            s0=z-c
            y=1.0_dp-s*s0*s0/(2.0_dp*z*b)
         else
            s5=s-1.0_dp
            s6=-2.0_dp*bc+b*s
            if(z<bc-s6/s5)then
               s7=z-bc
               y=-(s6/2.0_dp+s5*s5/s6*s7*(s7/2.0_dp+s6/s5))/z
            else
               y=0.0_dp
            end if
         end if
      end if
   end function lqq_weight

   elemental function hard_rejection_weight(x,cutoff) result(y)
      real(dp),intent(in)::x,cutoff
      real(dp)::y
      if(abs(x)<=cutoff) then
         y=1.0_dp
      else
         y=0.0_dp
      end if
   end function hard_rejection_weight

   elemental function smooth_weight(x,cutoff) result(y)
      real(dp),intent(in)::x,cutoff
      real(dp)::y,u
      u=abs(x)/cutoff
      if(u>=1.0_dp) then
         y=0.0_dp
      else
         y=(1.0_dp-u*u)**2
      end if
   end function smooth_weight
end module robustbase_psi
