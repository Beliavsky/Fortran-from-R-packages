! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from Risk 1.0 by Saralees Nadarajah and Stephen Chan.
! Copyright (c) 2017 Saralees Nadarajah and Stephen Chan.
module risk_measures
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use risk_kinds, only : dp
   use risk_math, only : quiet_nan
   use risk_distributions, only : continuous_distribution
   use risk_numerics, only : numeric_function, integrate, solve_bisection, &
      is_negative_infinity, is_positive_infinity
   implicit none
   private

   integer, parameter :: mode_first_moment = 1
   integer, parameter :: mode_second_moment = 2
   integer, parameter :: mode_survival = 3
   integer, parameter :: mode_cdf = 4
   integer, parameter :: mode_expectile_left = 5
   integer, parameter :: mode_expectile_right = 6
   integer, parameter :: mode_downside_square = 7
   integer, parameter :: mode_downside_power = 8
   integer, parameter :: mode_wang_upper = 9
   integer, parameter :: mode_wang_lower = 10
   integer, parameter :: mode_stone = 11
   integer, parameter :: mode_entropy = 12
   integer, parameter :: mode_density_power = 13
   integer, parameter :: mode_log_x = 14
   integer, parameter :: mode_x_power = 15
   integer, parameter :: mode_exp_cx = 16
   integer, parameter :: mode_pdf = 17
   integer, parameter :: mode_log_abs = 18
   integer, parameter :: mode_log_abs_square = 19
   integer, parameter :: mode_abs_power = 20
   integer, parameter :: mode_abs_double_power = 21
   integer, parameter :: mode_mean_abs_deviation = 22

   type, extends(numeric_function) :: distribution_integrand
      class(continuous_distribution), pointer :: dist => null()
      integer :: mode = 0
      real(dp) :: p1 = 0.0_dp
      real(dp) :: p2 = 0.0_dp
   contains
      procedure :: evaluate => evaluate_distribution_integrand
   end type distribution_integrand

   type, extends(numeric_function) :: expectile_balance
      class(continuous_distribution), pointer :: dist => null()
      real(dp) :: alpha = 0.5_dp
      real(dp) :: a = 0.0_dp
      real(dp) :: b = 0.0_dp
   contains
      procedure :: evaluate => evaluate_expectile_balance
   end type expectile_balance

   public :: varg, esg, tcm, expp, bvar, epsg, expect, expvar
   public :: omegag, sortinog, kappag, wangg1, wangg2
   public :: stoneg1, stoneg2, luceg1, luceg2, luceg3, luceg4
   public :: saring1, saring2, saring3
   public :: bkg1, bkg2, bkg3, bkg4

   interface varg
      module procedure varg_scalar, varg_vector
   end interface varg
   interface esg
      module procedure esg_scalar, esg_vector
   end interface esg
   interface tcm
      module procedure tcm_scalar, tcm_vector
   end interface tcm
   interface expp
      module procedure expp_scalar, expp_vector
   end interface expp
   interface bvar
      module procedure bvar_scalar, bvar_vector
   end interface bvar
   interface epsg
      module procedure epsg_scalar, epsg_vector
   end interface epsg
   interface expvar
      module procedure expvar_scalar, expvar_vector
   end interface expvar
   interface omegag
      module procedure omegag_scalar, omegag_vector
   end interface omegag
   interface sortinog
      module procedure sortinog_scalar, sortinog_vector
   end interface sortinog
   interface kappag
      module procedure kappag_scalar, kappag_vector
   end interface kappag
   interface wangg1
      module procedure wangg1_scalar, wangg1_vector
   end interface wangg1
   interface wangg2
      module procedure wangg2_scalar, wangg2_vector
   end interface wangg2
   interface bkg1
      module procedure bkg1_scalar, bkg1_vector
   end interface bkg1
   interface bkg2
      module procedure bkg2_scalar, bkg2_vector
   end interface bkg2
   interface bkg3
      module procedure bkg3_scalar, bkg3_vector
   end interface bkg3
   interface bkg4
      module procedure bkg4_scalar, bkg4_vector
   end interface bkg4

contains

   function evaluate_distribution_integrand(self, x) result(y)
      class(distribution_integrand), intent(in) :: self
      real(dp), intent(in) :: x
      real(dp) :: y
      real(dp) :: density, probability, exponent_value

      if (.not. associated(self%dist)) then
         y = quiet_nan()
         return
      end if

      select case (self%mode)
      case (mode_first_moment)
         y = x*self%dist%pdf(x)
      case (mode_second_moment)
         y = x*x*self%dist%pdf(x)
      case (mode_survival)
         y = 1.0_dp-self%dist%cdf(x)
      case (mode_cdf)
         y = self%dist%cdf(x)
      case (mode_expectile_left)
         y = max(self%p1-x,0.0_dp)**2*self%dist%pdf(x)
      case (mode_expectile_right)
         y = max(x-self%p1,0.0_dp)**2*self%dist%pdf(x)
      case (mode_downside_square)
         y = max(self%p1-x,0.0_dp)**2*self%dist%pdf(x)
      case (mode_downside_power)
         y = max(self%p1-x,0.0_dp)**self%p2*self%dist%pdf(x)
      case (mode_wang_upper)
         probability = max(0.0_dp,min(1.0_dp,1.0_dp-self%dist%cdf(x)))
         y = probability**self%p1-probability
      case (mode_wang_lower)
         probability = max(0.0_dp,min(1.0_dp,self%dist%cdf(x)))
         y = probability**self%p1-probability
      case (mode_stone)
         y = abs(x-self%p1)**self%p2*self%dist%pdf(x)
      case (mode_entropy)
         density = self%dist%pdf(x)
         if (density <= 0.0_dp) then
            y = 0.0_dp
         else
            y = log(density)*density
         end if
      case (mode_density_power)
         density = self%dist%pdf(x)
         exponent_value = 1.0_dp-self%p1
         if (density > 0.0_dp) then
            y = density**exponent_value
         else if (exponent_value > 0.0_dp) then
            y = 0.0_dp
         else if (abs(exponent_value) <= epsilon(1.0_dp)) then
            y = 1.0_dp
         else
            y = quiet_nan()
         end if
      case (mode_log_x)
         if (abs(x) <= tiny(1.0_dp)) then
            y = 0.0_dp
         else if (x < 0.0_dp) then
            y = quiet_nan()
         else
            y = log(x)*self%dist%pdf(x)
         end if
      case (mode_x_power)
         if (abs(x) <= tiny(1.0_dp) .and. self%p1 < 0.0_dp) then
            y = quiet_nan()
         else if (x < 0.0_dp .and. abs(self%p1-anint(self%p1)) > &
                  100.0_dp*epsilon(1.0_dp)) then
            y = quiet_nan()
         else
            y = x**self%p1*self%dist%pdf(x)
         end if
      case (mode_exp_cx)
         y = exp(self%p1*x)*self%dist%pdf(x)
      case (mode_pdf)
         y = self%dist%pdf(x)
      case (mode_log_abs)
         if (abs(x) <= tiny(1.0_dp)) then
            y = 0.0_dp
         else
            y = log(abs(x))*self%dist%pdf(x)
         end if
      case (mode_log_abs_square)
         if (abs(x) <= tiny(1.0_dp)) then
            y = 0.0_dp
         else
            y = log(abs(x))**2*self%dist%pdf(x)
         end if
      case (mode_abs_power)
         y = abs(x)**self%p1*self%dist%pdf(x)
      case (mode_abs_double_power)
         y = abs(x)**(2.0_dp*self%p1)*self%dist%pdf(x)
      case (mode_mean_abs_deviation)
         y = abs(x-self%p1)*self%dist%pdf(x)
      case default
         y = quiet_nan()
      end select
   end function evaluate_distribution_integrand

   function evaluate_expectile_balance(self, x) result(y)
      class(expectile_balance), intent(in) :: self
      real(dp), intent(in) :: x
      real(dp) :: y, left, right

      if (.not. associated(self%dist)) then
         y = quiet_nan()
         return
      end if
      left = integrate_mode(self%dist,mode_expectile_left,self%a,min(x,self%b),x,0.0_dp, &
                            1.0e-8_dp,1.0e-8_dp)
      right = integrate_mode(self%dist,mode_expectile_right,max(x,self%a),self%b,x,0.0_dp, &
                             1.0e-8_dp,1.0e-8_dp)
      y = self%alpha*left-(1.0_dp-self%alpha)*right
   end function evaluate_expectile_balance

   function integrate_mode(dist, mode, a, b, p1, p2, abs_tol, rel_tol) result(value)
      class(continuous_distribution), target, intent(in) :: dist
      integer, intent(in) :: mode
      real(dp), intent(in) :: a, b, p1, p2
      real(dp), intent(in), optional :: abs_tol, rel_tol
      real(dp) :: value
      type(distribution_integrand) :: fun
      real(dp) :: atol, rtol

      fun%dist => dist
      fun%mode = mode
      fun%p1 = p1
      fun%p2 = p2
      atol = 1.0e-9_dp
      rtol = 1.0e-9_dp
      if (present(abs_tol)) atol = abs_tol
      if (present(rel_tol)) rtol = rel_tol
      value = integrate(fun,a,b,atol,rtol)
   end function integrate_mode

   function split_at_zero_mode(dist, mode, a, b, p1, p2) result(value)
      class(continuous_distribution), target, intent(in) :: dist
      integer, intent(in) :: mode
      real(dp), intent(in) :: a, b, p1, p2
      real(dp) :: value

      if (a < 0.0_dp .and. b > 0.0_dp) then
         value = integrate_mode(dist,mode,a,0.0_dp,p1,p2)+ &
                 integrate_mode(dist,mode,0.0_dp,b,p1,p2)
      else
         value = integrate_mode(dist,mode,a,b,p1,p2)
      end if
   end function split_at_zero_mode

   function varg_scalar(dist, alpha) result(value)
      class(continuous_distribution), intent(in) :: dist
      real(dp), intent(in) :: alpha
      real(dp) :: value
      if (alpha < 0.0_dp .or. alpha > 1.0_dp) then
         value = quiet_nan()
      else
         value = dist%quantile(alpha)
      end if
   end function varg_scalar

   function varg_vector(dist, alpha) result(values)
      class(continuous_distribution), intent(in) :: dist
      real(dp), intent(in) :: alpha(:)
      real(dp) :: values(size(alpha))
      integer :: i
      do i = 1, size(alpha)
         values(i) = varg_scalar(dist,alpha(i))
      end do
   end function varg_vector

   function esg_scalar(dist, alpha) result(value)
      class(continuous_distribution), target, intent(in) :: dist
      real(dp), intent(in) :: alpha
      real(dp) :: value, q
      if (alpha <= 0.0_dp .or. alpha > 1.0_dp) then
         value = quiet_nan()
         return
      end if
      q = dist%quantile(alpha)
      value = integrate_mode(dist,mode_first_moment,dist%lower_bound(),q,0.0_dp,0.0_dp)/alpha
   end function esg_scalar

   function esg_vector(dist, alpha) result(values)
      class(continuous_distribution), intent(in) :: dist
      real(dp), intent(in) :: alpha(:)
      real(dp) :: values(size(alpha))
      integer :: i
      do i = 1, size(alpha)
         values(i) = esg_scalar(dist,alpha(i))
      end do
   end function esg_vector

   function tcm_scalar(dist, alpha) result(value)
      class(continuous_distribution), intent(in) :: dist
      real(dp), intent(in) :: alpha
      real(dp) :: value
      if (alpha < 0.0_dp .or. alpha >= 1.0_dp) then
         value = quiet_nan()
      else
         value = dist%quantile(0.5_dp*(1.0_dp+alpha))
      end if
   end function tcm_scalar

   function tcm_vector(dist, alpha) result(values)
      class(continuous_distribution), intent(in) :: dist
      real(dp), intent(in) :: alpha(:)
      real(dp) :: values(size(alpha))
      integer :: i
      do i = 1, size(alpha)
         values(i) = tcm_scalar(dist,alpha(i))
      end do
   end function tcm_vector

   function expp_scalar(dist, alpha, a, b) result(value)
      class(continuous_distribution), target, intent(in) :: dist
      real(dp), intent(in) :: alpha, a, b
      real(dp) :: value, lo, hi, flo, fhi, span
      integer :: i
      type(expectile_balance) :: balance

      if (alpha <= 0.0_dp .or. alpha >= 1.0_dp .or. a >= b) then
         value = quiet_nan()
         return
      end if
      if (is_negative_infinity(a)) then
         lo = dist%quantile(1.0e-8_dp)
      else
         lo = a
      end if
      if (is_positive_infinity(b)) then
         hi = dist%quantile(1.0_dp-1.0e-8_dp)
      else
         hi = b
      end if
      if (.not. ieee_is_finite(lo) .or. .not. ieee_is_finite(hi) .or. lo >= hi) then
         value = quiet_nan()
         return
      end if

      balance%dist => dist
      balance%alpha = alpha
      balance%a = a
      balance%b = b
      flo = balance%evaluate(lo)
      fhi = balance%evaluate(hi)
      span = max(1.0_dp,hi-lo)
      do i = 1, 20
         if (flo*fhi <= 0.0_dp) exit
         lo = lo-span
         hi = hi+span
         span = 2.0_dp*span
         flo = balance%evaluate(lo)
         fhi = balance%evaluate(hi)
      end do
      if (.not. ieee_is_finite(flo) .or. .not. ieee_is_finite(fhi) .or. flo*fhi > 0.0_dp) then
         value = quiet_nan()
      else
         value = solve_bisection(balance,lo,hi,1.0e-8_dp,1.0e-9_dp,160)
      end if
   end function expp_scalar

   function expp_vector(dist, alpha, a, b) result(values)
      class(continuous_distribution), intent(in) :: dist
      real(dp), intent(in) :: alpha(:)
      real(dp), intent(in) :: a, b
      real(dp) :: values(size(alpha))
      integer :: i
      do i = 1, size(alpha)
         values(i) = expp_scalar(dist,alpha(i),a,b)
      end do
   end function expp_vector

   function bvar_scalar(dist, alpha, a) result(value)
      class(continuous_distribution), target, intent(in) :: dist
      real(dp), intent(in) :: alpha, a
      real(dp) :: value, q
      if (alpha <= 0.0_dp .or. alpha > 1.0_dp) then
         value = quiet_nan()
         return
      end if
      q = dist%quantile(alpha)
      value = integrate_mode(dist,mode_first_moment,a,q,0.0_dp,0.0_dp)/alpha
   end function bvar_scalar

   function bvar_vector(dist, alpha, a) result(values)
      class(continuous_distribution), intent(in) :: dist
      real(dp), intent(in) :: alpha(:)
      real(dp), intent(in) :: a
      real(dp) :: values(size(alpha))
      integer :: i
      do i = 1, size(alpha)
         values(i) = bvar_scalar(dist,alpha(i),a)
      end do
   end function bvar_vector

   function epsg_scalar(dist, alpha) result(value)
      class(continuous_distribution), intent(in) :: dist
      real(dp), intent(in) :: alpha
      real(dp) :: value, q
      q = varg_scalar(dist,alpha)
      if (.not. ieee_is_finite(q) .or. abs(q) <= tiny(1.0_dp)) then
         value = quiet_nan()
      else
         value = (1.0_dp-alpha)*esg_scalar(dist,alpha)/q
      end if
   end function epsg_scalar

   function epsg_vector(dist, alpha) result(values)
      class(continuous_distribution), intent(in) :: dist
      real(dp), intent(in) :: alpha(:)
      real(dp) :: values(size(alpha))
      integer :: i
      do i = 1, size(alpha)
         values(i) = epsg_scalar(dist,alpha(i))
      end do
   end function epsg_vector

   function expect(dist, a, b) result(value)
      class(continuous_distribution), target, intent(in) :: dist
      real(dp), intent(in) :: a, b
      real(dp) :: value
      value = integrate_mode(dist,mode_first_moment,a,b,0.0_dp,0.0_dp)
   end function expect

   function expvar_scalar(dist, alpha, a, b) result(value)
      class(continuous_distribution), target, intent(in) :: dist
      real(dp), intent(in) :: alpha, a, b
      real(dp) :: value, mean_value, second_moment, variance
      mean_value = expect(dist,a,b)
      second_moment = integrate_mode(dist,mode_second_moment,a,b,0.0_dp,0.0_dp)
      variance = max(0.0_dp,second_moment-mean_value*mean_value)
      value = mean_value+alpha*sqrt(variance)
   end function expvar_scalar

   function expvar_vector(dist, alpha, a, b) result(values)
      class(continuous_distribution), intent(in) :: dist
      real(dp), intent(in) :: alpha(:)
      real(dp), intent(in) :: a, b
      real(dp) :: values(size(alpha))
      integer :: i
      do i = 1, size(alpha)
         values(i) = expvar_scalar(dist,alpha(i),a,b)
      end do
   end function expvar_vector

   function omegag_scalar(dist, alpha, a, b) result(value)
      class(continuous_distribution), target, intent(in) :: dist
      real(dp), intent(in) :: alpha, a, b
      real(dp) :: value, numerator, denominator
      numerator = integrate_mode(dist,mode_survival,alpha,b,0.0_dp,0.0_dp)
      denominator = integrate_mode(dist,mode_cdf,a,alpha,0.0_dp,0.0_dp)
      if (abs(denominator) <= tiny(1.0_dp)) then
         value = quiet_nan()
      else
         value = numerator/denominator
      end if
   end function omegag_scalar

   function omegag_vector(dist, alpha, a, b) result(values)
      class(continuous_distribution), intent(in) :: dist
      real(dp), intent(in) :: alpha(:)
      real(dp), intent(in) :: a, b
      real(dp) :: values(size(alpha))
      integer :: i
      do i = 1, size(alpha)
         values(i) = omegag_scalar(dist,alpha(i),a,b)
      end do
   end function omegag_vector

   function sortinog_scalar(dist, alpha, a, b) result(value)
      class(continuous_distribution), target, intent(in) :: dist
      real(dp), intent(in) :: alpha, a, b
      real(dp) :: value, mean_value, downside
      mean_value = expect(dist,a,b)
      downside = integrate_mode(dist,mode_downside_square,a,min(alpha,b),alpha,0.0_dp)
      if (downside <= 0.0_dp) then
         value = quiet_nan()
      else
         value = (mean_value-alpha)/sqrt(downside)
      end if
   end function sortinog_scalar

   function sortinog_vector(dist, alpha, a, b) result(values)
      class(continuous_distribution), intent(in) :: dist
      real(dp), intent(in) :: alpha(:)
      real(dp), intent(in) :: a, b
      real(dp) :: values(size(alpha))
      integer :: i
      do i = 1, size(alpha)
         values(i) = sortinog_scalar(dist,alpha(i),a,b)
      end do
   end function sortinog_vector

   function kappag_scalar(dist, alpha, n, a, b) result(value)
      class(continuous_distribution), target, intent(in) :: dist
      real(dp), intent(in) :: alpha, n, a, b
      real(dp) :: value, mean_value, lower_moment
      if (n <= 0.0_dp) then
         value = quiet_nan()
         return
      end if
      mean_value = expect(dist,a,b)
      lower_moment = integrate_mode(dist,mode_downside_power,a,min(alpha,b),alpha,n)
      if (lower_moment <= 0.0_dp) then
         value = quiet_nan()
      else
         value = (mean_value-alpha)*lower_moment**(-1.0_dp/n)
      end if
   end function kappag_scalar

   function kappag_vector(dist, alpha, n, a, b) result(values)
      class(continuous_distribution), intent(in) :: dist
      real(dp), intent(in) :: alpha(:)
      real(dp), intent(in) :: n, a, b
      real(dp) :: values(size(alpha))
      integer :: i
      do i = 1, size(alpha)
         values(i) = kappag_scalar(dist,alpha(i),n,a,b)
      end do
   end function kappag_vector

   function wangg1_scalar(dist, alpha, a, b) result(value)
      class(continuous_distribution), target, intent(in) :: dist
      real(dp), intent(in) :: alpha, a, b
      real(dp) :: value
      if (alpha <= 0.0_dp) then
         value = quiet_nan()
      else
         value = 0.5_dp*integrate_mode(dist,mode_wang_upper,a,b,alpha,0.0_dp)+ &
                 0.5_dp*integrate_mode(dist,mode_wang_lower,a,b,alpha,0.0_dp)
      end if
   end function wangg1_scalar

   function wangg1_vector(dist, alpha, a, b) result(values)
      class(continuous_distribution), intent(in) :: dist
      real(dp), intent(in) :: alpha(:)
      real(dp), intent(in) :: a, b
      real(dp) :: values(size(alpha))
      integer :: i
      do i = 1, size(alpha)
         values(i) = wangg1_scalar(dist,alpha(i),a,b)
      end do
   end function wangg1_vector

   function wangg2_scalar(dist, alpha, a, b) result(value)
      class(continuous_distribution), target, intent(in) :: dist
      real(dp), intent(in) :: alpha, a, b
      real(dp) :: value
      if (alpha <= 0.0_dp) then
         value = quiet_nan()
      else
         value = integrate_mode(dist,mode_wang_upper,a,b,alpha,0.0_dp)
      end if
   end function wangg2_scalar

   function wangg2_vector(dist, alpha, a, b) result(values)
      class(continuous_distribution), intent(in) :: dist
      real(dp), intent(in) :: alpha(:)
      real(dp), intent(in) :: a, b
      real(dp) :: values(size(alpha))
      integer :: i
      do i = 1, size(alpha)
         values(i) = wangg2_scalar(dist,alpha(i),a,b)
      end do
   end function wangg2_vector

   function stoneg1(dist, x0, k, a, b) result(value)
      class(continuous_distribution), target, intent(in) :: dist
      real(dp), intent(in) :: x0, k, a, b
      real(dp) :: value
      value = integrate_mode(dist,mode_stone,a,b,x0,k)
   end function stoneg1

   function stoneg2(dist, x0, k, a, b) result(value)
      class(continuous_distribution), intent(in) :: dist
      real(dp), intent(in) :: x0, k, a, b
      real(dp) :: value
      if (k <= 0.0_dp) then
         value = quiet_nan()
      else
         value = stoneg1(dist,x0,k,a,b)**(1.0_dp/k)
      end if
   end function stoneg2

   function luceg1(dist, a, b, aa, bb) result(value)
      class(continuous_distribution), target, intent(in) :: dist
      real(dp), intent(in) :: a, b, aa, bb
      real(dp) :: value
      value = bb-aa*integrate_mode(dist,mode_entropy,a,b,0.0_dp,0.0_dp)
   end function luceg1

   function luceg2(dist, a, b, aa, bb) result(value)
      class(continuous_distribution), target, intent(in) :: dist
      real(dp), intent(in) :: a, b, aa, bb
      real(dp) :: value
      value = aa*integrate_mode(dist,mode_density_power,a,b,bb,0.0_dp)
   end function luceg2

   function luceg3(dist, a, b, aa, bb) result(value)
      class(continuous_distribution), target, intent(in) :: dist
      real(dp), intent(in) :: a, b, aa, bb
      real(dp) :: value
      value = bb+aa*integrate_mode(dist,mode_log_x,a,b,0.0_dp,0.0_dp)
   end function luceg3

   function luceg4(dist, a, b, aa, bb) result(value)
      class(continuous_distribution), target, intent(in) :: dist
      real(dp), intent(in) :: a, b, aa, bb
      real(dp) :: value
      value = aa*integrate_mode(dist,mode_x_power,a,b,bb,0.0_dp)
   end function luceg4

   function saring1(dist, a, b, k, c) result(value)
      class(continuous_distribution), target, intent(in) :: dist
      real(dp), intent(in) :: a, b, k, c
      real(dp) :: value
      value = k*integrate_mode(dist,mode_exp_cx,a,b,c,0.0_dp)
   end function saring1

   function saring2(dist, a, b, aa, bb1, bb2) result(value)
      class(continuous_distribution), target, intent(in) :: dist
      real(dp), intent(in) :: a, b, aa, bb1, bb2
      real(dp) :: value, positive_mass, negative_mass, first_log, second_log
      positive_mass = 0.0_dp
      negative_mass = 0.0_dp
      if (b > 0.0_dp) positive_mass = integrate_mode(dist,mode_pdf,max(a,0.0_dp),b,0.0_dp,0.0_dp)
      if (a < 0.0_dp) negative_mass = integrate_mode(dist,mode_pdf,a,min(b,0.0_dp),0.0_dp,0.0_dp)
      first_log = split_at_zero_mode(dist,mode_log_abs,a,b,0.0_dp,0.0_dp)
      second_log = split_at_zero_mode(dist,mode_log_abs_square,a,b,0.0_dp,0.0_dp)
      value = bb1*positive_mass+bb2*negative_mass+aa*first_log- &
              0.5_dp*aa*aa*(second_log-first_log*first_log)
   end function saring2

   function saring3(dist, a, b, aa, bb1, bb2) result(value)
      class(continuous_distribution), target, intent(in) :: dist
      real(dp), intent(in) :: a, b, aa, bb1, bb2
      real(dp) :: value, positive_part, negative_part, first_abs, second_abs
      if (abs(aa-1.0_dp) <= 100.0_dp*epsilon(1.0_dp)) then
         value = quiet_nan()
         return
      end if
      positive_part = 0.0_dp
      negative_part = 0.0_dp
      if (b > 0.0_dp) positive_part = integrate_mode(dist,mode_abs_power,max(a,0.0_dp),b,aa,0.0_dp)
      if (a < 0.0_dp) negative_part = integrate_mode(dist,mode_abs_power,a,min(b,0.0_dp),aa,0.0_dp)
      first_abs = integrate_mode(dist,mode_abs_power,a,b,aa,0.0_dp)
      second_abs = integrate_mode(dist,mode_abs_double_power,a,b,aa,0.0_dp)
      value = bb1*positive_part+bb2*negative_part+ &
              aa*second_abs/(2.0_dp*(aa-1.0_dp))- &
              0.5_dp*(second_abs-first_abs*first_abs)
   end function saring3

   function bkg1_scalar(dist, alpha, a, b) result(value)
      class(continuous_distribution), intent(in) :: dist
      real(dp), intent(in) :: alpha, a, b
      real(dp) :: value
      value = varg_scalar(dist,alpha)-expect(dist,a,b)
   end function bkg1_scalar

   function bkg1_vector(dist, alpha, a, b) result(values)
      class(continuous_distribution), intent(in) :: dist
      real(dp), intent(in) :: alpha(:)
      real(dp), intent(in) :: a, b
      real(dp) :: values(size(alpha))
      integer :: i
      do i = 1, size(alpha)
         values(i) = bkg1_scalar(dist,alpha(i),a,b)
      end do
   end function bkg1_vector

   function bkg2_scalar(dist, alpha, a, b) result(value)
      class(continuous_distribution), intent(in) :: dist
      real(dp), intent(in) :: alpha, a, b
      real(dp) :: value
      value = esg_scalar(dist,alpha)-expect(dist,a,b)
   end function bkg2_scalar

   function bkg2_vector(dist, alpha, a, b) result(values)
      class(continuous_distribution), intent(in) :: dist
      real(dp), intent(in) :: alpha(:)
      real(dp), intent(in) :: a, b
      real(dp) :: values(size(alpha))
      integer :: i
      do i = 1, size(alpha)
         values(i) = bkg2_scalar(dist,alpha(i),a,b)
      end do
   end function bkg2_vector

   function bkg3_scalar(dist, alpha, a, b, beta) result(value)
      class(continuous_distribution), target, intent(in) :: dist
      real(dp), intent(in) :: alpha, a, b, beta
      real(dp) :: value, mean_value, mean_absolute_deviation
      mean_value = expect(dist,a,b)
      mean_absolute_deviation = integrate_mode(dist,mode_mean_abs_deviation,a,b,mean_value,0.0_dp)
      value = esg_scalar(dist,alpha)-beta*mean_absolute_deviation
   end function bkg3_scalar

   function bkg3_vector(dist, alpha, a, b, beta) result(values)
      class(continuous_distribution), intent(in) :: dist
      real(dp), intent(in) :: alpha(:)
      real(dp), intent(in) :: a, b, beta
      real(dp) :: values(size(alpha))
      integer :: i
      do i = 1, size(alpha)
         values(i) = bkg3_scalar(dist,alpha(i),a,b,beta)
      end do
   end function bkg3_vector

   function bkg4_scalar(dist, alpha, a, b, beta) result(value)
      class(continuous_distribution), target, intent(in) :: dist
      real(dp), intent(in) :: alpha, a, b, beta
      real(dp) :: value, mean_value, mean_absolute_deviation
      mean_value = expect(dist,a,b)
      mean_absolute_deviation = integrate_mode(dist,mode_mean_abs_deviation,a,b,mean_value,0.0_dp)
      value = varg_scalar(dist,alpha)+esg_scalar(dist,alpha)-beta*mean_absolute_deviation
   end function bkg4_scalar

   function bkg4_vector(dist, alpha, a, b, beta) result(values)
      class(continuous_distribution), intent(in) :: dist
      real(dp), intent(in) :: alpha(:)
      real(dp), intent(in) :: a, b, beta
      real(dp) :: values(size(alpha))
      integer :: i
      do i = 1, size(alpha)
         values(i) = bkg4_scalar(dist,alpha(i),a,b,beta)
      end do
   end function bkg4_vector

end module risk_measures
