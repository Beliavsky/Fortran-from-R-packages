! SPDX-License-Identifier: MIT
module gkwdist_ad
   use gkwdist_kinds, only : dp
   use gkwdist_math, only : digamma_fn, trigamma_fn, log1mexp, expm1_stable
   implicit none
   private
   integer, parameter, public :: ad_dim = 5

   type, public :: ad2
      real(dp) :: v = 0.0_dp
      real(dp) :: g(ad_dim) = 0.0_dp
      real(dp) :: h(ad_dim,ad_dim) = 0.0_dp
   end type ad2

   public :: ad_const, ad_var, ad_log_gamma, ad_log1mexp
   public :: operator(+), operator(-), operator(*), operator(/), log, exp

   interface operator(+)
      module procedure add_aa, add_ar, add_ra
   end interface
   interface operator(-)
      module procedure sub_aa, sub_ar, sub_ra, neg_a
   end interface
   interface operator(*)
      module procedure mul_aa, mul_ar, mul_ra
   end interface
   interface operator(/)
      module procedure div_aa, div_ar, div_ra
   end interface
   interface log
      module procedure log_a
   end interface
   interface exp
      module procedure exp_a
   end interface

contains

   pure function ad_const(x) result(a)
      real(dp), intent(in) :: x
      type(ad2) :: a
      a%v = x
   end function ad_const

   pure function ad_var(x, idx) result(a)
      real(dp), intent(in) :: x
      integer, intent(in) :: idx
      type(ad2) :: a
      a%v = x
      if (idx >= 1 .and. idx <= ad_dim) a%g(idx) = 1.0_dp
   end function ad_var

   pure function add_aa(a,b) result(c)
      type(ad2), intent(in) :: a,b
      type(ad2) :: c
      c%v = a%v+b%v; c%g=a%g+b%g; c%h=a%h+b%h
   end function add_aa
   pure function add_ar(a,b) result(c)
      type(ad2), intent(in) :: a; real(dp), intent(in) :: b
      type(ad2) :: c
      c=a; c%v=a%v+b
   end function add_ar
   pure function add_ra(a,b) result(c)
      real(dp), intent(in) :: a; type(ad2), intent(in) :: b
      type(ad2) :: c
      c=b; c%v=a+b%v
   end function add_ra

   pure function sub_aa(a,b) result(c)
      type(ad2), intent(in) :: a,b
      type(ad2) :: c
      c%v=a%v-b%v; c%g=a%g-b%g; c%h=a%h-b%h
   end function sub_aa
   pure function sub_ar(a,b) result(c)
      type(ad2), intent(in) :: a; real(dp), intent(in) :: b
      type(ad2) :: c
      c=a; c%v=a%v-b
   end function sub_ar
   pure function sub_ra(a,b) result(c)
      real(dp), intent(in) :: a; type(ad2), intent(in) :: b
      type(ad2) :: c
      c%v=a-b%v; c%g=-b%g; c%h=-b%h
   end function sub_ra
   pure function neg_a(a) result(c)
      type(ad2), intent(in) :: a
      type(ad2) :: c
      c%v=-a%v; c%g=-a%g; c%h=-a%h
   end function neg_a

   pure function mul_aa(a,b) result(c)
      type(ad2), intent(in) :: a,b
      type(ad2) :: c
      integer :: i,j
      c%v=a%v*b%v
      c%g=a%g*b%v + b%g*a%v
      c%h=a%h*b%v + b%h*a%v
      do i=1,ad_dim
         do j=1,ad_dim
            c%h(i,j)=c%h(i,j)+a%g(i)*b%g(j)+b%g(i)*a%g(j)
         end do
      end do
   end function mul_aa
   pure function mul_ar(a,b) result(c)
      type(ad2), intent(in) :: a; real(dp), intent(in) :: b
      type(ad2) :: c
      c%v=a%v*b; c%g=a%g*b; c%h=a%h*b
   end function mul_ar
   pure function mul_ra(a,b) result(c)
      real(dp), intent(in) :: a; type(ad2), intent(in) :: b
      type(ad2) :: c
      c%v=a*b%v; c%g=a*b%g; c%h=a*b%h
   end function mul_ra

   pure function inv_a(a) result(c)
      type(ad2), intent(in) :: a
      type(ad2) :: c
      real(dp) :: f1, f2
      integer :: i,j
      c%v=1.0_dp/a%v
      f1=-1.0_dp/(a%v*a%v)
      f2=2.0_dp/(a%v*a%v*a%v)
      c%g=f1*a%g
      c%h=f1*a%h
      do i=1,ad_dim
         do j=1,ad_dim
            c%h(i,j)=c%h(i,j)+f2*a%g(i)*a%g(j)
         end do
      end do
   end function inv_a
   pure function div_aa(a,b) result(c)
      type(ad2), intent(in) :: a,b
      type(ad2) :: c
      c=a*inv_a(b)
   end function div_aa
   pure function div_ar(a,b) result(c)
      type(ad2), intent(in) :: a; real(dp), intent(in) :: b
      type(ad2) :: c
      c=a*(1.0_dp/b)
   end function div_ar
   pure function div_ra(a,b) result(c)
      real(dp), intent(in) :: a; type(ad2), intent(in) :: b
      type(ad2) :: c
      c=a*inv_a(b)
   end function div_ra

   pure function log_a(a) result(c)
      type(ad2), intent(in) :: a
      type(ad2) :: c
      real(dp) :: f1, f2
      integer :: i,j
      c%v=log(a%v)
      f1=1.0_dp/a%v; f2=-1.0_dp/(a%v*a%v)
      c%g=f1*a%g; c%h=f1*a%h
      do i=1,ad_dim
         do j=1,ad_dim
            c%h(i,j)=c%h(i,j)+f2*a%g(i)*a%g(j)
         end do
      end do
   end function log_a

   pure function exp_a(a) result(c)
      type(ad2), intent(in) :: a
      type(ad2) :: c
      integer :: i,j
      c%v=exp(a%v)
      c%g=c%v*a%g; c%h=c%v*a%h
      do i=1,ad_dim
         do j=1,ad_dim
            c%h(i,j)=c%h(i,j)+c%v*a%g(i)*a%g(j)
         end do
      end do
   end function exp_a

   pure function ad_log1mexp(a) result(c)
      type(ad2), intent(in) :: a
      type(ad2) :: c
      real(dp) :: t, denom, f1, f2
      integer :: i,j
      c%v=log1mexp(a%v)
      t=exp(a%v)
      denom=-expm1_stable(a%v)
      f1=-t/denom
      f2=-t/(denom*denom)
      c%g=f1*a%g; c%h=f1*a%h
      do i=1,ad_dim
         do j=1,ad_dim
            c%h(i,j)=c%h(i,j)+f2*a%g(i)*a%g(j)
         end do
      end do
   end function ad_log1mexp

   pure function ad_log_gamma(a) result(c)
      type(ad2), intent(in) :: a
      type(ad2) :: c
      real(dp) :: f1, f2
      integer :: i,j
      c%v=log_gamma(a%v)
      f1=digamma_fn(a%v)
      f2=trigamma_fn(a%v)
      c%g=f1*a%g; c%h=f1*a%h
      do i=1,ad_dim
         do j=1,ad_dim
            c%h(i,j)=c%h(i,j)+f2*a%g(i)*a%g(j)
         end do
      end do
   end function ad_log_gamma

end module gkwdist_ad
