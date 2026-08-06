! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from ghyp 1.6.5 by Marc Weibel, David Luethi, and Henriette-Elise Breymann.
module ghyp_special
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
   use ghyp_kinds, only : dp, pi
   implicit none
   private

   public :: normal_pdf, normal_cdf, normal_quantile
   public :: regularized_beta, student_cdf
   public :: digamma_fn, log_bessel_k, bessel_k
   public :: gauss_legendre_integral, gauss_legendre_rule
   public :: log_sum_exp

   real(dp), allocatable, save :: nodes128(:), weights128(:)
   real(dp), allocatable, save :: nodes192(:), weights192(:)
   real(dp), allocatable, save :: nodes256(:), weights256(:)

   abstract interface
      function scalar_function(x) result(y)
         import dp
         real(dp), intent(in) :: x
         real(dp) :: y
      end function scalar_function
   end interface

contains

   pure elemental function normal_pdf(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value
      value = exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
   end function normal_pdf

   pure elemental function normal_cdf(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value
      value = 0.5_dp*erfc(-x/sqrt(2.0_dp))
   end function normal_cdf

   pure function normal_quantile(p) result(x)
      real(dp), intent(in) :: p
      real(dp) :: x
      real(dp) :: q, r
      real(dp), parameter :: a(6) = [ &
         -3.969683028665376e1_dp, 2.209460984245205e2_dp, &
         -2.759285104469687e2_dp, 1.383577518672690e2_dp, &
         -3.066479806614716e1_dp, 2.506628277459239_dp]
      real(dp), parameter :: b(5) = [ &
         -5.447609879822406e1_dp, 1.615858368580409e2_dp, &
         -1.556989798598866e2_dp, 6.680131188771972e1_dp, &
         -1.328068155288572e1_dp]
      real(dp), parameter :: c(6) = [ &
         -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
         -2.400758277161838_dp, -2.549732539343734_dp, &
          4.374664141464968_dp, 2.938163982698783_dp]
      real(dp), parameter :: d(4) = [ &
          7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
          2.445134137142996_dp, 3.754408661907416_dp]
      real(dp), parameter :: plow = 0.02425_dp
      real(dp), parameter :: phigh = 1.0_dp - plow

      if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else if (p < plow) then
         q = sqrt(-2.0_dp*log(p))
         x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
             ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      else if (p <= phigh) then
         q = p - 0.5_dp
         r = q*q
         x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/ &
             (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
      else
         q = sqrt(-2.0_dp*log(1.0_dp-p))
         x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
              ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      end if
      if (abs(x) < huge(1.0_dp)/2.0_dp) then
         x = x - (normal_cdf(x)-p)/max(normal_pdf(x), tiny(1.0_dp))
      end if
   end function normal_quantile

   pure function betacf(a, b, x) result(value)
      real(dp), intent(in) :: a, b, x
      real(dp) :: value
      integer :: m, m2
      real(dp) :: aa, c, d, del, h, qab, qam, qap
      real(dp), parameter :: eps = 3.0e-14_dp
      real(dp), parameter :: fpmin = 1.0e-300_dp

      qab = a+b
      qap = a+1.0_dp
      qam = a-1.0_dp
      c = 1.0_dp
      d = 1.0_dp-qab*x/qap
      if (abs(d) < fpmin) d = fpmin
      d = 1.0_dp/d
      h = d
      do m = 1, 300
         m2 = 2*m
         aa = real(m,dp)*(b-real(m,dp))*x/((qam+real(m2,dp))*(a+real(m2,dp)))
         d = 1.0_dp+aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp+aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         h = h*d*c
         aa = -(a+real(m,dp))*(qab+real(m,dp))*x/ &
              ((a+real(m2,dp))*(qap+real(m2,dp)))
         d = 1.0_dp+aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp+aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         del = d*c
         h = h*del
         if (abs(del-1.0_dp) <= eps) exit
      end do
      value = h
   end function betacf

   pure function regularized_beta(x, a, b) result(value)
      real(dp), intent(in) :: x, a, b
      real(dp) :: value, bt
      if (x <= 0.0_dp) then
         value = 0.0_dp
         return
      else if (x >= 1.0_dp) then
         value = 1.0_dp
         return
      else if (a <= 0.0_dp .or. b <= 0.0_dp) then
         value = ieee_value(1.0_dp, ieee_quiet_nan)
         return
      end if
      bt = exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log(1.0_dp-x))
      if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
         value = bt*betacf(a,b,x)/a
      else
         value = 1.0_dp-bt*betacf(b,a,1.0_dp-x)/b
      end if
      value = min(1.0_dp,max(0.0_dp,value))
   end function regularized_beta

   pure function student_cdf(x, nu) result(value)
      real(dp), intent(in) :: x, nu
      real(dp) :: value, z
      if (nu <= 0.0_dp) then
         value = ieee_value(1.0_dp, ieee_quiet_nan)
         return
      end if
      if (abs(x) <= tiny(1.0_dp)) then
         value = 0.5_dp
         return
      end if
      z = nu/(nu+x*x)
      if (x > 0.0_dp) then
         value = 1.0_dp-0.5_dp*regularized_beta(z,0.5_dp*nu,0.5_dp)
      else
         value = 0.5_dp*regularized_beta(z,0.5_dp*nu,0.5_dp)
      end if
   end function student_cdf

   pure function digamma_fn(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value, y, r, f
      if (x <= 0.0_dp) then
         value = ieee_value(1.0_dp, ieee_quiet_nan)
         return
      end if
      y = x
      r = 0.0_dp
      do while (y < 8.0_dp)
         r = r-1.0_dp/y
         y = y+1.0_dp
      end do
      f = 1.0_dp/(y*y)
      value = r+log(y)-0.5_dp/y-f*(1.0_dp/12.0_dp- &
              f*(1.0_dp/120.0_dp-f*(1.0_dp/252.0_dp-f/240.0_dp)))
   end function digamma_fn

   subroutine compute_gl_rule(nn, nodes, weights)
      integer, intent(in) :: nn
      real(dp), allocatable, intent(out) :: nodes(:), weights(:)
      integer :: i, j, m
      real(dp) :: z, z1, p1, p2, p3, pp
      allocate(nodes(nn),weights(nn))
      m = (nn+1)/2
      do i = 1, m
         z = cos(pi*(real(i,dp)-0.25_dp)/(real(nn,dp)+0.5_dp))
         do
            p1 = 1.0_dp
            p2 = 0.0_dp
            do j = 1, nn
               p3 = p2
               p2 = p1
               p1 = ((2.0_dp*real(j,dp)-1.0_dp)*z*p2- &
                    (real(j,dp)-1.0_dp)*p3)/real(j,dp)
            end do
            pp = real(nn,dp)*(z*p1-p2)/(z*z-1.0_dp)
            z1 = z
            z = z1-p1/pp
            if (abs(z-z1) <= 4.0_dp*epsilon(z)) exit
         end do
         nodes(i) = -z
         nodes(nn+1-i) = z
         weights(i) = 2.0_dp/((1.0_dp-z*z)*pp*pp)
         weights(nn+1-i) = weights(i)
      end do
   end subroutine compute_gl_rule

   subroutine gauss_legendre_rule(nn, nodes, weights)
      integer, intent(in) :: nn
      real(dp), allocatable, intent(out) :: nodes(:), weights(:)
      select case(nn)
      case(128)
         if (.not. allocated(nodes128)) call compute_gl_rule(128,nodes128,weights128)
         allocate(nodes(128),weights(128))
         nodes=nodes128
         weights=weights128
      case(192)
         if (.not. allocated(nodes192)) call compute_gl_rule(192,nodes192,weights192)
         allocate(nodes(192),weights(192))
         nodes=nodes192
         weights=weights192
      case(256)
         if (.not. allocated(nodes256)) call compute_gl_rule(256,nodes256,weights256)
         allocate(nodes(256),weights(256))
         nodes=nodes256
         weights=weights256
      case default
         call compute_gl_rule(nn,nodes,weights)
      end select
   end subroutine gauss_legendre_rule

   recursive function gauss_legendre_integral(f, a, b, n) result(value)
      procedure(scalar_function) :: f
      real(dp), intent(in) :: a, b
      integer, intent(in), optional :: n
      real(dp) :: value
      real(dp), allocatable :: local_nodes(:), local_weights(:)
      real(dp) :: xm, xl
      integer :: i, nn
      nn = 128
      if (present(n)) nn = max(8,n)
      select case(nn)
      case(128)
         if (.not. allocated(nodes128)) call compute_gl_rule(128,nodes128,weights128)
         xm=0.5_dp*(a+b)
         xl=0.5_dp*(b-a)
         value=0.0_dp
         do i=1,128
            value=value+weights128(i)*f(xm+xl*nodes128(i))
         end do
      case(192)
         if (.not. allocated(nodes192)) call compute_gl_rule(192,nodes192,weights192)
         xm=0.5_dp*(a+b)
         xl=0.5_dp*(b-a)
         value=0.0_dp
         do i=1,192
            value=value+weights192(i)*f(xm+xl*nodes192(i))
         end do
      case(256)
         if (.not. allocated(nodes256)) call compute_gl_rule(256,nodes256,weights256)
         xm=0.5_dp*(a+b)
         xl=0.5_dp*(b-a)
         value=0.0_dp
         do i=1,256
            value=value+weights256(i)*f(xm+xl*nodes256(i))
         end do
      case default
         call compute_gl_rule(nn,local_nodes,local_weights)
         xm=0.5_dp*(a+b)
         xl=0.5_dp*(b-a)
         value=0.0_dp
         do i=1,nn
            value=value+local_weights(i)*f(xm+xl*local_nodes(i))
         end do
      end select
      value=xl*value
   end function gauss_legendre_integral

   pure function log_cosh(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value, ax
      ax = abs(x)
      value = ax+log(1.0_dp+exp(-2.0_dp*ax))-log(2.0_dp)
   end function log_cosh

   function log_bessel_k(nu, x) result(value)
      real(dp), intent(in) :: nu, x
      real(dp) :: value, integral, tmax, t, ly, half
      real(dp), allocatable :: nodes(:), weights(:)
      integer :: i
      if (x <= 0.0_dp) then
         value = huge(1.0_dp)
         return
      end if
      if (abs(abs(nu)-0.5_dp) <= 8.0_dp*epsilon(1.0_dp)) then
         value = 0.5_dp*log(pi/(2.0_dp*x))-x
         return
      end if
      tmax = max(10.0_dp, log(2.0_dp*(55.0_dp+10.0_dp*abs(nu))/x+2.0_dp))
      tmax = min(tmax,30.0_dp)
      call gauss_legendre_rule(192,nodes,weights)
      half=0.5_dp*tmax
      integral=0.0_dp
      do i=1,size(nodes)
         t=half*(nodes(i)+1.0_dp)
         ly=-x*(cosh(t)-1.0_dp)+log_cosh(nu*t)
         if(ly>log(tiny(1.0_dp)).and.ly<log(huge(1.0_dp))) &
            integral=integral+weights(i)*exp(ly)
      end do
      integral=half*integral
      if (integral <= 0.0_dp .or. .not. ieee_is_finite(integral)) then
         value = 0.5_dp*log(pi/(2.0_dp*x))-x + &
            log(max(tiny(1.0_dp),1.0_dp+(4.0_dp*nu*nu-1.0_dp)/(8.0_dp*x)))
      else
         value = log(integral)-x
      end if
   end function log_bessel_k

   function bessel_k(nu, x) result(value)
      real(dp), intent(in) :: nu, x
      real(dp) :: value, lv
      lv = log_bessel_k(nu,x)
      if (lv > log(huge(1.0_dp))) then
         value = huge(1.0_dp)
      else if (lv < log(tiny(1.0_dp))) then
         value = 0.0_dp
      else
         value = exp(lv)
      end if
   end function bessel_k

   pure function log_sum_exp(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value, m
      m = maxval(x)
      if (.not. ieee_is_finite(m)) then
         value = m
      else
         value = m+log(sum(exp(x-m)))
      end if
   end function log_sum_exp

end module ghyp_special
