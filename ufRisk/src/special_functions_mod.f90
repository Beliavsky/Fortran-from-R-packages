! SPDX-License-Identifier: GPL-3.0-only
module special_functions_mod
   use kind_mod, only: dp
   implicit none
   private
   public :: regularized_gamma_q, regularized_beta
contains
   pure real(dp) function regularized_gamma_q(a, x) result(q)
      real(dp), intent(in) :: a, x
      integer, parameter :: maxit=400
      real(dp), parameter :: eps=2.0e-14_dp, fpmin=1.0e-300_dp
      real(dp) :: sumv, del, ap, b, c, d, h, an
      integer :: i
      if (a <= 0.0_dp .or. x < 0.0_dp) then
         q = 0.0_dp; return
      else if (x <= tiny(1.0_dp)) then
         q = 1.0_dp; return
      end if
      if (x < a+1.0_dp) then
         ap=a; sumv=1.0_dp/a; del=sumv
         do i=1,maxit
            ap=ap+1.0_dp; del=del*x/ap; sumv=sumv+del
            if (abs(del) <= abs(sumv)*eps) exit
         end do
         q = 1.0_dp-sumv*exp(-x+a*log(x)-log_gamma(a))
      else
         b=x+1.0_dp-a; c=1.0_dp/fpmin; d=1.0_dp/max(b,fpmin); h=d
         do i=1,maxit
            an=-real(i,dp)*(real(i,dp)-a)
            b=b+2.0_dp; d=an*d+b; if(abs(d)<fpmin)d=fpmin
            c=b+an/c; if(abs(c)<fpmin)c=fpmin
            d=1.0_dp/d; del=d*c; h=h*del
            if(abs(del-1.0_dp)<eps)exit
         end do
         q=exp(-x+a*log(x)-log_gamma(a))*h
      end if
      q=max(0.0_dp,min(1.0_dp,q))
   end function regularized_gamma_q

   pure real(dp) function regularized_beta(x,a,b) result(value)
      real(dp), intent(in) :: x,a,b
      real(dp) :: bt
      if (a<=0.0_dp .or. b<=0.0_dp) then
         value=0.0_dp; return
      else if (x<=0.0_dp) then
         value=0.0_dp; return
      else if (x>=1.0_dp) then
         value=1.0_dp; return
      end if
      bt=exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log(1.0_dp-x))
      if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
         value=bt*beta_cf(x,a,b)/a
      else
         value=1.0_dp-bt*beta_cf(1.0_dp-x,b,a)/b
      end if
      value=max(0.0_dp,min(1.0_dp,value))
   end function regularized_beta

   pure real(dp) function beta_cf(x,a,b) result(h)
      real(dp), intent(in) :: x,a,b
      integer, parameter :: maxit=400
      real(dp), parameter :: eps=2.0e-14_dp, fpmin=1.0e-300_dp
      real(dp) :: qab,qap,qam,c,d,aa,del
      integer :: m,m2
      qab=a+b; qap=a+1.0_dp; qam=a-1.0_dp
      c=1.0_dp; d=1.0_dp-qab*x/qap; if(abs(d)<fpmin)d=fpmin
      d=1.0_dp/d; h=d
      do m=1,maxit
         m2=2*m
         aa=real(m,dp)*(b-real(m,dp))*x/((qam+real(m2,dp))*(a+real(m2,dp)))
         d=1.0_dp+aa*d; if(abs(d)<fpmin)d=fpmin
         c=1.0_dp+aa/c; if(abs(c)<fpmin)c=fpmin
         d=1.0_dp/d; h=h*d*c
         aa=-(a+real(m,dp))*(qab+real(m,dp))*x/ &
            ((a+real(m2,dp))*(qap+real(m2,dp)))
         d=1.0_dp+aa*d; if(abs(d)<fpmin)d=fpmin
         c=1.0_dp+aa/c; if(abs(c)<fpmin)c=fpmin
         d=1.0_dp/d; del=d*c; h=h*del
         if(abs(del-1.0_dp)<eps)exit
      end do
   end function beta_cf
end module special_functions_mod
