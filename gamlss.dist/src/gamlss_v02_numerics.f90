! Numerical helpers for extended gamlss.dist families.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_v02_numerics
   use gamlss_kinds, only : dp, pi
   implicit none
   private
   public :: adaptive_integral, bisection_root, bessel_k, log_bessel_k

   abstract interface
      function scalar_fun(x) result(y)
         import dp
         real(dp), intent(in) :: x
         real(dp) :: y
      end function scalar_fun
   end interface
contains

   recursive function adaptive_integral(f,a,b,tol,max_depth) result(val)
      procedure(scalar_fun) :: f
      real(dp), intent(in) :: a,b
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: max_depth
      real(dp) :: val,fa,fb,fc,c,s,eps
      integer :: depth
      if (a==b) then
         val=0.0_dp
         return
      end if
      eps=1.0e-9_dp
      if (present(tol)) eps=max(tol,1.0e-14_dp)
      depth=20
      if (present(max_depth)) depth=max(4,max_depth)
      c=0.5_dp*(a+b)
      fa=f(a)
      fb=f(b)
      fc=f(c)
      s=(b-a)*(fa+4.0_dp*fc+fb)/6.0_dp
      val=asimpson(f,a,b,fa,fb,fc,s,eps,depth)
   end function adaptive_integral

   recursive function asimpson(f,a,b,fa,fb,fc,s,eps,depth) result(v)
      procedure(scalar_fun) :: f
      real(dp), intent(in) :: a,b,fa,fb,fc,s,eps
      integer, intent(in) :: depth
      real(dp) :: v,c,d,e,fd,fe,s1,s2
      c=0.5_dp*(a+b)
      d=0.5_dp*(a+c)
      e=0.5_dp*(c+b)
      fd=f(d)
      fe=f(e)
      s1=(c-a)*(fa+4.0_dp*fd+fc)/6.0_dp
      s2=(b-c)*(fc+4.0_dp*fe+fb)/6.0_dp
      if (depth<=0 .or. abs(s1+s2-s)<=15.0_dp*eps) then
         v=s1+s2+(s1+s2-s)/15.0_dp
      else
         v=asimpson(f,a,c,fa,fc,fd,s1,0.5_dp*eps,depth-1) &
           +asimpson(f,c,b,fc,fb,fe,s2,0.5_dp*eps,depth-1)
      end if
   end function asimpson

   function bisection_root(f,target,lo_in,hi_in,tol,max_iter) result(root)
      procedure(scalar_fun) :: f
      real(dp), intent(in) :: target,lo_in,hi_in
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: max_iter
      real(dp) :: root,lo,hi,mid,fm,eps
      integer :: i,nit
      lo=lo_in
      hi=hi_in
      eps=1.0e-9_dp
      if (present(tol)) eps=max(tol,1.0e-14_dp)
      nit=120
      if (present(max_iter)) nit=max(10,max_iter)
      do i=1,nit
         mid=0.5_dp*(lo+hi)
         fm=f(mid)
         if (fm<target) then
            lo=mid
         else
            hi=mid
         end if
         if (abs(hi-lo)<=eps*max(1.0_dp,abs(mid))) exit
      end do
      root=0.5_dp*(lo+hi)
   end function bisection_root

   function log_bessel_k(x,nu) result(lk)
      real(dp), intent(in) :: x,nu
      real(dp) :: lk,av,lead,c1,c2,c3,c4,c5,kval,mu4
      if (x<=0.0_dp) then
         lk=huge(1.0_dp)
         return
      end if
      av=abs(nu)
      if (x>35.0_dp) then
         lead=0.5_dp*log(pi/(2.0_dp*x))-x
         mu4=4.0_dp*av*av
         c1=(mu4-1.0_dp)/(8.0_dp*x)
         c2=(mu4-1.0_dp)*(mu4-9.0_dp)/(128.0_dp*x*x)
         c3=(mu4-1.0_dp)*(mu4-9.0_dp)*(mu4-25.0_dp)/(3072.0_dp*x*x*x)
         c4=(mu4-1.0_dp)*(mu4-9.0_dp)*(mu4-25.0_dp)*(mu4-49.0_dp) &
            /(98304.0_dp*x**4)
         c5=(mu4-1.0_dp)*(mu4-9.0_dp)*(mu4-25.0_dp)*(mu4-49.0_dp)*(mu4-81.0_dp) &
            /(3932160.0_dp*x**5)
         lk=lead+log(max(1.0_dp+c1+c2+c3+c4+c5,tiny(1.0_dp)))
         return
      end if
      if (x<1.0e-5_dp .and. av>1.0e-6_dp) then
         lk=log(0.5_dp)+log_gamma(av)-av*log(0.5_dp*x)
         return
      end if
      kval=bessel_k_integral(x,av)
      if (kval<=0.0_dp) then
         lk=-huge(1.0_dp)
      else
         lk=log(kval)
      end if
   end function log_bessel_k

   function bessel_k(x,nu) result(k)
      real(dp), intent(in) :: x,nu
      real(dp) :: k,lk
      lk=log_bessel_k(x,nu)
      if (lk>log(huge(1.0_dp))) then
         k=huge(1.0_dp)
      else if (lk<log(tiny(1.0_dp))) then
         k=0.0_dp
      else
         k=exp(lk)
      end if
   end function bessel_k

   function bessel_k_integral(x,nu) result(k)
      real(dp), intent(in) :: x,nu
      real(dp) :: k,tmax
      tmax=2.0_dp
      do while (x*cosh(tmax)-nu*tmax<38.0_dp .and. tmax<30.0_dp)
         tmax=tmax+1.0_dp
      end do
      k=adaptive_integral(integrand,0.0_dp,tmax,2.0e-11_dp,24)
   contains
      function integrand(t) result(v)
         real(dp), intent(in) :: t
         real(dp) :: v,a,b,m
         a=-x*cosh(t)+nu*t
         b=-x*cosh(t)-nu*t
         m=max(a,b)
         if (m<-745.0_dp) then
            v=0.0_dp
         else
            v=0.5_dp*exp(m)*(exp(a-m)+exp(b-m))
         end if
      end function integrand
   end function bessel_k_integral

end module gamlss_v02_numerics
