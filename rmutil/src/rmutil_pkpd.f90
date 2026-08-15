! rmutil computational translation
! Copyright (C) 1998-2001 J.K. Lindsey
! Copyright (C) 2026 OpenAI (modern Fortran translation)
! SPDX-License-Identifier: GPL-2.0-or-later
module rmutil_pkpd
   use rmutil_kinds, only : dp
   implicit none
   private
   public :: mu1_0o1c, mu1_1o1c, mu1_1o2c, mu1_1o2cl, mu1_1o2cc
   public :: mu2_0o1c, mu2_0o2c1, mu2_0o2c2, mu2_1o1c
   public :: mu2_0o1cfp, mu2_0o2c1fp, mu2_0o2c2fp, mu2_1o1cfp
contains

   function mu1_0o1c(p,times,dose,end_time) result(y)
      real(dp),intent(in)::p(:),times(:)
      real(dp),intent(in),optional::dose,end_time
      real(dp),allocatable::y(:)
      real(dp)::ke,d,e
      integer::i
      d=1.0_dp; if(present(dose))d=dose
      e=0.5_dp; if(present(end_time))e=end_time
      ke=exp(p(2)); allocate(y(size(times)))
      do i=1,size(times)
         if(times(i)<=e)then
            y(i)=d/(exp(p(1))*ke)*(1.0_dp-exp(-ke*times(i)))
         else
            y(i)=d/(exp(p(1))*ke)*(1.0_dp-exp(-ke*e))*exp(-ke*(times(i)-e))
         end if
      end do
   end function mu1_0o1c

   function mu1_1o1c(p,times,dose) result(y)
      real(dp),intent(in)::p(:),times(:); real(dp),intent(in),optional::dose
      real(dp),allocatable::y(:); real(dp)::ka,ke,d
      d=1.0_dp; if(present(dose))d=dose
      ka=exp(p(2)); ke=exp(p(3)); allocate(y(size(times)))
      y=exp(p(2)-p(1))*d/(ka-ke)*(exp(-ke*times)-exp(-ka*times))
   end function mu1_1o1c

   function mu1_1o2c(p,times,dose) result(y)
      real(dp),intent(in)::p(:),times(:); real(dp),intent(in),optional::dose
      real(dp),allocatable::y(:); real(dp)::ka,ke,k12,d
      d=1.0_dp; if(present(dose))d=dose
      ka=exp(p(2)); ke=exp(p(3)); k12=exp(p(4)); allocate(y(size(times)))
      y=ka*k12*exp(-p(1))*d/(k12-ka)*((exp(-ka*times)-exp(-ke*times))/(ke-ka) - &
         (exp(-k12*times)-exp(-ke*times))/(ke-k12))
   end function mu1_1o2c

   function mu1_1o2cl(p,times,dose) result(y)
      real(dp),intent(in)::p(:),times(:); real(dp),intent(in),optional::dose
      real(dp),allocatable::y(:); real(dp)::ka,ke,d
      d=1.0_dp; if(present(dose))d=dose
      ka=exp(p(2)); ke=exp(p(3)); allocate(y(size(times)))
      y=ka*ka*exp(-p(1))*d/(ka-ke)*((exp(-ka*times)-exp(-ke*times))/(ke-ka) - &
         times*exp(-ka*times))
   end function mu1_1o2cl

   function mu1_1o2cc(p,times,dose) result(y)
      real(dp),intent(in)::p(:),times(:); real(dp),intent(in),optional::dose
      real(dp),allocatable::y(:)
      real(dp)::ka,ke,k12,k21,beta,alpha,d
      d=1.0_dp; if(present(dose))d=dose
      ka=exp(p(2)); ke=exp(p(3)); k12=exp(p(4)); k21=exp(p(5))
      beta=0.5_dp*(k12+k21+ke-sqrt((k12+k21+ke)**2-4.0_dp*k21*ke))
      alpha=(k21*ke)/beta; allocate(y(size(times)))
      y=exp(p(2)-p(1))*d*((k21-alpha)*exp(-alpha*times)/((ka-alpha)*(beta-alpha)) + &
         (k21-beta)*exp(-beta*times)/((ka-beta)*(alpha-beta)) + &
         (k21-ka)*exp(-ka*times)/((beta-ka)*(beta-ka)))
   end function mu1_1o2cc

   function mu2_0o1c(p,times,ind,dose,end_time) result(y)
      real(dp),intent(in)::p(:),times(:),ind(:)
      real(dp),intent(in),optional::dose,end_time
      real(dp),allocatable::y(:)
      real(dp)::vp,kpm,kp,kem,vm,kemp,d,e,t,tmp1,tmp2,g1,g2,cend,cexp
      real(dp)::cmend,tmp3,tmp4
      integer::i
      if(size(ind)/=size(times))error stop "mu2_0o1c: ind/times size mismatch"
      d=1.0_dp; if(present(dose))d=dose; e=0.5_dp; if(present(end_time))e=end_time
      vp=exp(p(1)); kpm=exp(p(3)); kp=exp(p(2))+kpm; kem=exp(p(4)); vm=exp(p(5))
      kemp=kem-kp; g1=exp(-kp*e); g2=exp(-kem*e); cend=(1.0_dp-g1)/(vp*kp)
      cmend=kpm/(kp*kem*vm)*(1.0_dp+kp/kemp*g2-kem/kemp*g1)
      tmp3=cend*kpm*vp/(kemp*vm)
      allocate(y(size(times)))
      do i=1,size(times)
         t=times(i); tmp1=exp(-kp*t); tmp2=kpm/(kp*kem*vm)
         if(t<=e)then
            y(i)=d*(ind(i)*(1.0_dp-tmp1)/(vp*kp) + (1.0_dp-ind(i))*tmp2* &
               (1.0_dp+kp/kemp*exp(-kem*t)-kem/kemp*tmp1))
         else
            cexp=exp(-kp*(t-e)); tmp4=g2/g1*cexp*tmp3 + &
               (cmend-tmp3)*exp(-kem*(t-e))/g2
            y(i)=d*(ind(i)*cend*cexp+(1.0_dp-ind(i))*tmp4)
         end if
      end do
   end function mu2_0o1c

   function mu2_0o2c1(p,times,ind,dose,end_time) result(y)
      real(dp),intent(in)::p(:),times(:),ind(:)
      real(dp),intent(in),optional::dose,end_time
      real(dp),allocatable::y(:)
      real(dp)::vp,kp12,kp21,kpm,kp,kem,d,e,t,tmp1,l1,l2
      real(dp)::t10,t13,t2,t3,t4,t5,t6,t7,t8,t9,t11,t12,t14,parent,metab
      integer::i
      d=1.0_dp; if(present(dose))d=dose; e=0.5_dp; if(present(end_time))e=end_time
      vp=exp(p(1)); kp12=exp(p(3)); kp21=exp(p(4)); kpm=exp(p(5)); kp=exp(p(2))+kpm
      kem=exp(p(6)); tmp1=sqrt((kp+kp12+kp21)**2-4.0_dp*kp21*kp)
      l1=0.5_dp*(kp+kp12+kp21+tmp1); l2=0.5_dp*(kp+kp12+kp21-tmp1)
      t13=exp(-kem*e); t3=(1.0_dp-t13)/kem; t4=l1-kp21; t5=l2-kp21
      t8=exp(-l1*e); t9=exp(-l2*e); allocate(y(size(times)))
      do i=1,size(times)
         t=times(i); t10=exp(-kem*t); t2=(1.0_dp-t10)/kem
         t6=exp(-l1*t); t7=exp(-l2*t); t11=t6/t8; t12=t7/t9; t14=t10/t13
         if(t<=e)then
            parent=t4*(1.0_dp-t6)/l1-t5*(1.0_dp-t7)/l2
            metab=kpm*(t4*(t2-(t6-t10)/(kem-l1))/l1 - &
               t5*(t2-(t7-t10)/(kem-l2))/l2)
         else
            parent=t4*(1.0_dp-t8)*t11/l1-t5*(1.0_dp-t9)*t12/l2
            metab=kpm*((t4*(t3-(t8-t13)/(kem-l1))/l1 - &
               t5*(t3-(t9-t13)/(kem-l2))/l2)*t14 + &
               t4*(1.0_dp-t8)*(t14-t11)/(l1*(l1-kem)) - &
               t5*(1.0_dp-t9)*(t14-t12)/(l2*(l2-kem)))
         end if
         y(i)=d/(vp*tmp1)*(ind(i)*parent+(1.0_dp-ind(i))*metab)
      end do
   end function mu2_0o2c1

   function mu2_0o2c2(p,times,ind,dose,end_time) result(y)
      real(dp),intent(in)::p(:),times(:),ind(:)
      real(dp),intent(in),optional::dose,end_time
      real(dp),allocatable::y(:)
      real(dp)::vp,kp12,kp21,kpm,kp,kem,km12,km21,d,e,t
      real(dp)::a,l1,l2,b,r1,r2,t2,t3,t6,t7,t11,t12,t13,t14
      real(dp)::t17,t18,t19,t20,t21,t22,t23,t24,parent,metab
      real(dp)::e1,e2,e3,e4,c1,c2,c3,c4
      integer::i
      d=1.0_dp; if(present(dose))d=dose; e=0.5_dp; if(present(end_time))e=end_time
      vp=exp(p(1)); kp12=exp(p(3)); kp21=exp(p(4)); kpm=exp(p(5)); kp=exp(p(2))+kpm
      kem=exp(p(6)); km12=exp(p(7)); km21=exp(p(8))
      a=sqrt((kp+kp12+kp21)**2-4.0_dp*kp21*kp)
      l1=0.5_dp*(kp+kp12+kp21+a); l2=0.5_dp*(kp+kp12+kp21-a)
      t2=l1-kp21; t3=l2-kp21; t6=exp(-l1*e); t7=exp(-l2*e)
      b=sqrt((kem+km12+km21)**2-4.0_dp*km21*kem)
      r1=0.5_dp*(kem+km12+km21+b); r2=0.5_dp*(kem+km12+km21-b)
      t11=kem-r1; t12=kem-r2; t13=l1-r2; t14=l1-r1
      t17=l1-km21; t18=l2-km21; t19=l2-r1; t20=l2-r2
      t21=exp(-r1*e); t22=exp(-r2*e); t23=km21-l1; t24=km21-l2
      allocate(y(size(times)))
      do i=1,size(times)
         t=times(i); e1=exp(-l1*t); e2=exp(-l2*t); e3=exp(-r1*t); e4=exp(-r2*t)
         if(t<=e)then
            parent=t2*(1.0_dp-e1)/(l1*a)-t3*(1.0_dp-e2)/(l2*a)
            c1=t2/a*((t11*e4/t13-t12*e3/t14)/(kem*b)+t17*e1/(l1*t14*t13)+1.0_dp/(kem*l1))
            c2=t3/a*((t12*e3/t19-t11*e4/t20)/(kem*(-b))+t18*e2/(l2*t19*t20)+1.0_dp/(kem*l2))
            metab=kpm*(c1-c2)
         else
            parent=t2*(1.0_dp-t6)*(e1/t6)/(l1*a)-t3*(1.0_dp-t7)*(e2/t7)/(l2*a)
            c1=(-t2/(a*b)*(t12*t21/(kem*t14)+r2/(kem*l1)+t23/(l1*t14)) + &
               t3/(a*b)*(t12*t21/(kem*t19)+r2/(kem*l2)+t24/(l2*t19)))*(e3/t21)
            c2=(t2/(a*b)*(t11*t22/(kem*t13)+r1/(kem*l1)+t23/(l1*t13)) - &
               t3/(a*b)*(t11*t22/(kem*t20)+r1/(kem*l2)+t24/(l2*t20)))*(e4/t22)
            c3=t2*t23*(1.0_dp-t6)*(e1/t6)/(l1*a*t14*t13)
            c4=t3*t24*(1.0_dp-t7)*(e2/t7)/(l2*a*t19*t20)
            metab=kpm*(c1+c2+c3-c4)
         end if
         y(i)=d/vp*(ind(i)*parent+(1.0_dp-ind(i))*metab)
      end do
   end function mu2_0o2c2

   function mu2_1o1c(p,times,ind,dose) result(y)
      real(dp),intent(in)::p(:),times(:),ind(:); real(dp),intent(in),optional::dose
      real(dp),allocatable::y(:); real(dp)::kap,kep,kem,d
      d=1.0_dp; if(present(dose))d=dose
      kap=exp(p(2)); kep=exp(p(3)); kem=exp(p(5)); allocate(y(size(times)))
      y=kap*exp(p(1))*d/(kap-kep)*(ind*(exp(-kep*times)-exp(-kap*times)) + &
         (1.0_dp-ind)*exp(p(4))*(exp(-kap*times)/(kap-kem)-exp(-kep*times)/(kep-kem) + &
         (1.0_dp/(kep-kem)-1.0_dp/(kap-kem))*exp(-kem*times)))
   end function mu2_1o1c

   function mu2_0o1cfp(p,times,ind,dose,end_time) result(y)
      real(dp),intent(in)::p(:),times(:),ind(:)
      real(dp),intent(in),optional::dose,end_time
      real(dp),allocatable::y(:)
      real(dp)::vp,kpm,kp,kem,vm,kemp,d,e,lpfp,t,tmp1,tmp2,tmp3,g1,g2,cend,cexp
      real(dp)::cmend,tmp4,parent,metab,direct
      integer::i
      d=1.0_dp; if(present(dose))d=dose; e=0.5_dp; if(present(end_time))e=end_time
      vp=exp(p(1)); kpm=exp(p(3)); kp=exp(p(2))+kpm; kem=exp(p(4)); vm=exp(p(5))
      kemp=kem-kp; tmp3=kpm/(kp*kem*vm); g1=exp(-kp*e); g2=exp(-kem*e)
      cend=(1.0_dp-g1)/(vp*kp); cmend=tmp3*(1.0_dp+kp/kemp*g2-kem/kemp*g1)
      tmp4=cend*kpm*vp/(kemp*vm); lpfp=exp(p(6)); lpfp=lpfp/(1.0_dp+lpfp)
      allocate(y(size(times)))
      do i=1,size(times)
         t=times(i); tmp1=exp(-kp*t); tmp2=exp(-kem*t)
         if(t<=e)then
            parent=(1.0_dp-tmp1)/(vp*kp)
            metab=tmp3*(1.0_dp+kp/kemp*tmp2-kem/kemp*tmp1)
         else
            cexp=exp(-kp*(t-e)); parent=cend*cexp
            metab=g2/g1*cexp*tmp4+(cmend-tmp4)*exp(-kem*(t-e))/g2
         end if
         if(t<=0.5_dp)then
            direct=(1.0_dp-tmp2)/(vm*kem)
         else
            direct=(1.0_dp-exp(-kem*0.5_dp))/(vm*kem)*exp(-kem*(t-0.5_dp))
         end if
         y(i)=d*(ind(i)*parent*lpfp+(1.0_dp-ind(i))*(direct*(1.0_dp-lpfp)+metab*lpfp))
      end do
   end function mu2_0o1cfp

   function mu2_0o2c1fp(p,times,ind,dose,end_time) result(y)
      real(dp),intent(in)::p(:),times(:),ind(:)
      real(dp),intent(in),optional::dose,end_time
      real(dp),allocatable::y(:),base(:)
      real(dp)::d,e,lpfp,vp,kem,direct
      integer::i
      d=1.0_dp; if(present(dose))d=dose; e=0.5_dp; if(present(end_time))e=end_time
      lpfp=exp(p(7)); lpfp=lpfp/(1.0_dp+lpfp); vp=exp(p(1)); kem=exp(p(6))
      base=mu2_0o2c1(p(1:6),times,ind,1.0_dp,e); allocate(y(size(times)))
      do i=1,size(times)
         if(times(i)<=e)then
            direct=(1.0_dp-exp(-kem*times(i)))/kem
         else
            direct=(1.0_dp-exp(-kem*e))/kem*exp(-kem*(times(i)-e))
         end if
         y(i)=d*(base(i)*lpfp + (1.0_dp-ind(i))*direct/vp*(1.0_dp-lpfp))
      end do
   end function mu2_0o2c1fp

   function mu2_0o2c2fp(p,times,ind,dose,end_time) result(y)
      real(dp),intent(in)::p(:),times(:),ind(:)
      real(dp),intent(in),optional::dose,end_time
      real(dp),allocatable::y(:),base(:)
      real(dp)::d,e,lpfp,vp,kem,direct
      integer::i
      d=1.0_dp; if(present(dose))d=dose; e=0.5_dp; if(present(end_time))e=end_time
      lpfp=exp(p(9)); lpfp=lpfp/(1.0_dp+lpfp); vp=exp(p(1)); kem=exp(p(6))
      base=mu2_0o2c2(p(1:8),times,ind,1.0_dp,e); allocate(y(size(times)))
      do i=1,size(times)
         if(times(i)<=e)then
            direct=(1.0_dp-exp(-kem*times(i)))/kem
         else
            direct=(1.0_dp-exp(-kem*e))/kem*exp(-kem*(times(i)-e))
         end if
         y(i)=d*(base(i)*lpfp+(1.0_dp-ind(i))*direct/vp*(1.0_dp-lpfp))
      end do
   end function mu2_0o2c2fp

   function mu2_1o1cfp(p,times,ind,dose) result(y)
      real(dp),intent(in)::p(:),times(:),ind(:); real(dp),intent(in),optional::dose
      real(dp),allocatable::y(:)
      real(dp)::kap,kf,kep,kmfp,kem,lpfp,d
      d=1.0_dp; if(present(dose))d=dose
      kap=exp(p(2)); kf=exp(p(4)); kep=exp(p(3))+kf; kmfp=exp(p(5)); kem=exp(p(6))
      lpfp=exp(p(7)); lpfp=lpfp/(1.0_dp+lpfp); allocate(y(size(times)))
      y=exp(-p(1))*d*(ind*(exp(-kep*times)-exp(-kap*times))/(kap-kep)*lpfp + &
         (1.0_dp-ind)*((exp(-kmfp*times)-exp(-kap*times))/(kap-kmfp)*(1.0_dp-lpfp) + &
         kf*(exp(-kap*times)/(kap-kem)-exp(-kep*times)/(kep-kem) + &
         (1.0_dp/(kep-kem)-1.0_dp/(kap-kem))*exp(-kem*times))/(kap-kep)*lpfp))
   end function mu2_1o1cfp

end module rmutil_pkpd
