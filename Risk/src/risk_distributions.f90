! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from Risk 1.0 by Saralees Nadarajah and Stephen Chan.
! Copyright (c) 2017 Saralees Nadarajah and Stephen Chan.
module risk_distributions
   use risk_kinds, only : dp
   use risk_math, only : pi, normal_pdf_std, normal_cdf_std, &
      normal_quantile_std, student_t_cdf_std, quiet_nan
   implicit none
   private

   type, abstract, public :: continuous_distribution
   contains
      procedure(pdf_method), deferred :: pdf
      procedure(cdf_method), deferred :: cdf
      procedure(quantile_method), deferred :: quantile
      procedure(bound_method), deferred :: lower_bound
      procedure(bound_method), deferred :: upper_bound
   end type continuous_distribution

   abstract interface
      function pdf_method(self, x) result(y)
         import :: continuous_distribution, dp
         class(continuous_distribution), intent(in) :: self
         real(dp), intent(in) :: x
         real(dp) :: y
      end function pdf_method

      function cdf_method(self, x) result(p)
         import :: continuous_distribution, dp
         class(continuous_distribution), intent(in) :: self
         real(dp), intent(in) :: x
         real(dp) :: p
      end function cdf_method

      function quantile_method(self, p) result(x)
         import :: continuous_distribution, dp
         class(continuous_distribution), intent(in) :: self
         real(dp), intent(in) :: p
         real(dp) :: x
      end function quantile_method

      function bound_method(self) result(x)
         import :: continuous_distribution, dp
         class(continuous_distribution), intent(in) :: self
         real(dp) :: x
      end function bound_method

      function scalar_callback(x) result(y)
         import :: dp
         real(dp), intent(in) :: x
         real(dp) :: y
      end function scalar_callback
   end interface

   type, extends(continuous_distribution), public :: normal_distribution
      real(dp) :: mu = 0.0_dp
      real(dp) :: sigma = 1.0_dp
   contains
      procedure :: pdf => normal_pdf
      procedure :: cdf => normal_cdf
      procedure :: quantile => normal_quantile
      procedure :: lower_bound => normal_lower
      procedure :: upper_bound => normal_upper
   end type normal_distribution

   type, extends(continuous_distribution), public :: lognormal_distribution
      real(dp) :: meanlog = 0.0_dp
      real(dp) :: sdlog = 1.0_dp
   contains
      procedure :: pdf => lognormal_pdf
      procedure :: cdf => lognormal_cdf
      procedure :: quantile => lognormal_quantile
      procedure :: lower_bound => lognormal_lower
      procedure :: upper_bound => lognormal_upper
   end type lognormal_distribution

   type, extends(continuous_distribution), public :: uniform_distribution
      real(dp) :: a = 0.0_dp
      real(dp) :: b = 1.0_dp
   contains
      procedure :: pdf => uniform_pdf
      procedure :: cdf => uniform_cdf
      procedure :: quantile => uniform_quantile
      procedure :: lower_bound => uniform_lower
      procedure :: upper_bound => uniform_upper
   end type uniform_distribution

   type, extends(continuous_distribution), public :: exponential_distribution
      real(dp) :: rate = 1.0_dp
   contains
      procedure :: pdf => exponential_pdf
      procedure :: cdf => exponential_cdf
      procedure :: quantile => exponential_quantile
      procedure :: lower_bound => exponential_lower
      procedure :: upper_bound => exponential_upper
   end type exponential_distribution

   type, extends(continuous_distribution), public :: logistic_distribution
      real(dp) :: mu = 0.0_dp
      real(dp) :: scale = 1.0_dp
   contains
      procedure :: pdf => logistic_pdf
      procedure :: cdf => logistic_cdf
      procedure :: quantile => logistic_quantile
      procedure :: lower_bound => logistic_lower
      procedure :: upper_bound => logistic_upper
   end type logistic_distribution

   type, extends(continuous_distribution), public :: student_t_distribution
      real(dp) :: nu = 5.0_dp
      real(dp) :: mu = 0.0_dp
      real(dp) :: scale = 1.0_dp
   contains
      procedure :: pdf => student_t_pdf
      procedure :: cdf => student_t_cdf
      procedure :: quantile => student_t_quantile
      procedure :: lower_bound => student_t_lower
      procedure :: upper_bound => student_t_upper
   end type student_t_distribution

   type, extends(continuous_distribution), public :: callback_distribution
      procedure(scalar_callback), pointer, nopass :: pdf_proc => null()
      procedure(scalar_callback), pointer, nopass :: cdf_proc => null()
      procedure(scalar_callback), pointer, nopass :: quantile_proc => null()
      real(dp) :: lower = -huge(1.0_dp)
      real(dp) :: upper = huge(1.0_dp)
   contains
      procedure :: pdf => callback_pdf
      procedure :: cdf => callback_cdf
      procedure :: quantile => callback_quantile
      procedure :: lower_bound => callback_lower
      procedure :: upper_bound => callback_upper
      procedure :: initialize => initialize_callback_distribution
   end type callback_distribution

contains

   function normal_lower(self) result(x)
      class(normal_distribution), intent(in) :: self
      real(dp) :: x
      x = -huge(1.0_dp)+0.0_dp*self%mu
   end function normal_lower

   function normal_upper(self) result(x)
      class(normal_distribution), intent(in) :: self
      real(dp) :: x
      x = huge(1.0_dp)+0.0_dp*self%mu
   end function normal_upper

   function lognormal_lower(self) result(x)
      class(lognormal_distribution), intent(in) :: self
      real(dp) :: x
      x = 0.0_dp*self%meanlog
   end function lognormal_lower

   function lognormal_upper(self) result(x)
      class(lognormal_distribution), intent(in) :: self
      real(dp) :: x
      x = huge(1.0_dp)+0.0_dp*self%meanlog
   end function lognormal_upper

   function exponential_lower(self) result(x)
      class(exponential_distribution), intent(in) :: self
      real(dp) :: x
      x = 0.0_dp*self%rate
   end function exponential_lower

   function exponential_upper(self) result(x)
      class(exponential_distribution), intent(in) :: self
      real(dp) :: x
      x = huge(1.0_dp)+0.0_dp*self%rate
   end function exponential_upper

   function logistic_lower(self) result(x)
      class(logistic_distribution), intent(in) :: self
      real(dp) :: x
      x = -huge(1.0_dp)+0.0_dp*self%mu
   end function logistic_lower

   function logistic_upper(self) result(x)
      class(logistic_distribution), intent(in) :: self
      real(dp) :: x
      x = huge(1.0_dp)+0.0_dp*self%mu
   end function logistic_upper

   function student_t_lower(self) result(x)
      class(student_t_distribution), intent(in) :: self
      real(dp) :: x
      x = -huge(1.0_dp)+0.0_dp*self%mu
   end function student_t_lower

   function student_t_upper(self) result(x)
      class(student_t_distribution), intent(in) :: self
      real(dp) :: x
      x = huge(1.0_dp)+0.0_dp*self%mu
   end function student_t_upper

   function normal_pdf(self, x) result(y)
      class(normal_distribution), intent(in) :: self
      real(dp), intent(in) :: x
      real(dp) :: y, z
      if (self%sigma <= 0.0_dp) then
         y = quiet_nan()
      else
         z = (x-self%mu)/self%sigma
         y = normal_pdf_std(z)/self%sigma
      end if
   end function normal_pdf

   function normal_cdf(self, x) result(p)
      class(normal_distribution), intent(in) :: self
      real(dp), intent(in) :: x
      real(dp) :: p
      if (self%sigma <= 0.0_dp) then
         p = quiet_nan()
      else
         p = normal_cdf_std((x-self%mu)/self%sigma)
      end if
   end function normal_cdf

   function normal_quantile(self, p) result(x)
      class(normal_distribution), intent(in) :: self
      real(dp), intent(in) :: p
      real(dp) :: x
      if (self%sigma <= 0.0_dp) then
         x = quiet_nan()
      else
         x = self%mu+self%sigma*normal_quantile_std(p)
      end if
   end function normal_quantile

   function lognormal_pdf(self, x) result(y)
      class(lognormal_distribution), intent(in) :: self
      real(dp), intent(in) :: x
      real(dp) :: y, z
      if (self%sdlog <= 0.0_dp) then
         y = quiet_nan()
      else if (x <= 0.0_dp) then
         y = 0.0_dp
      else
         z = (log(x)-self%meanlog)/self%sdlog
         y = normal_pdf_std(z)/(x*self%sdlog)
      end if
   end function lognormal_pdf

   function lognormal_cdf(self, x) result(p)
      class(lognormal_distribution), intent(in) :: self
      real(dp), intent(in) :: x
      real(dp) :: p
      if (self%sdlog <= 0.0_dp) then
         p = quiet_nan()
      else if (x <= 0.0_dp) then
         p = 0.0_dp
      else
         p = normal_cdf_std((log(x)-self%meanlog)/self%sdlog)
      end if
   end function lognormal_cdf

   function lognormal_quantile(self, p) result(x)
      class(lognormal_distribution), intent(in) :: self
      real(dp), intent(in) :: p
      real(dp) :: x
      if (self%sdlog <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
         x = quiet_nan()
      else if (p <= 0.0_dp) then
         x = 0.0_dp
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else
         x = exp(self%meanlog+self%sdlog*normal_quantile_std(p))
      end if
   end function lognormal_quantile

   function uniform_pdf(self, x) result(y)
      class(uniform_distribution), intent(in) :: self
      real(dp), intent(in) :: x
      real(dp) :: y
      if (self%b <= self%a) then
         y = quiet_nan()
      else if (x < self%a .or. x > self%b) then
         y = 0.0_dp
      else
         y = 1.0_dp/(self%b-self%a)
      end if
   end function uniform_pdf

   function uniform_cdf(self, x) result(p)
      class(uniform_distribution), intent(in) :: self
      real(dp), intent(in) :: x
      real(dp) :: p
      if (self%b <= self%a) then
         p = quiet_nan()
      else if (x <= self%a) then
         p = 0.0_dp
      else if (x >= self%b) then
         p = 1.0_dp
      else
         p = (x-self%a)/(self%b-self%a)
      end if
   end function uniform_cdf

   function uniform_quantile(self, p) result(x)
      class(uniform_distribution), intent(in) :: self
      real(dp), intent(in) :: p
      real(dp) :: x
      if (self%b <= self%a .or. p < 0.0_dp .or. p > 1.0_dp) then
         x = quiet_nan()
      else
         x = self%a+p*(self%b-self%a)
      end if
   end function uniform_quantile

   function uniform_lower(self) result(x)
      class(uniform_distribution), intent(in) :: self
      real(dp) :: x
      x = self%a
   end function uniform_lower

   function uniform_upper(self) result(x)
      class(uniform_distribution), intent(in) :: self
      real(dp) :: x
      x = self%b
   end function uniform_upper

   function exponential_pdf(self, x) result(y)
      class(exponential_distribution), intent(in) :: self
      real(dp), intent(in) :: x
      real(dp) :: y
      if (self%rate <= 0.0_dp) then
         y = quiet_nan()
      else if (x < 0.0_dp) then
         y = 0.0_dp
      else
         y = self%rate*exp(-self%rate*x)
      end if
   end function exponential_pdf

   function exponential_cdf(self, x) result(p)
      class(exponential_distribution), intent(in) :: self
      real(dp), intent(in) :: x
      real(dp) :: p
      if (self%rate <= 0.0_dp) then
         p = quiet_nan()
      else if (x <= 0.0_dp) then
         p = 0.0_dp
      else
         p = 1.0_dp-exp(-self%rate*x)
      end if
   end function exponential_cdf

   function exponential_quantile(self, p) result(x)
      class(exponential_distribution), intent(in) :: self
      real(dp), intent(in) :: p
      real(dp) :: x
      if (self%rate <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
         x = quiet_nan()
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else
         x = -log(1.0_dp-p)/self%rate
      end if
   end function exponential_quantile

   function logistic_pdf(self, x) result(y)
      class(logistic_distribution), intent(in) :: self
      real(dp), intent(in) :: x
      real(dp) :: y, z, e
      if (self%scale <= 0.0_dp) then
         y = quiet_nan()
      else
         z = (x-self%mu)/self%scale
         if (z >= 0.0_dp) then
            e = exp(-z)
         else
            e = exp(z)
         end if
         y = e/(self%scale*(1.0_dp+e)**2)
      end if
   end function logistic_pdf

   function logistic_cdf(self, x) result(p)
      class(logistic_distribution), intent(in) :: self
      real(dp), intent(in) :: x
      real(dp) :: p, z
      if (self%scale <= 0.0_dp) then
         p = quiet_nan()
      else
         z = (x-self%mu)/self%scale
         if (z >= 0.0_dp) then
            p = 1.0_dp/(1.0_dp+exp(-z))
         else
            p = exp(z)/(1.0_dp+exp(z))
         end if
      end if
   end function logistic_cdf

   function logistic_quantile(self, p) result(x)
      class(logistic_distribution), intent(in) :: self
      real(dp), intent(in) :: p
      real(dp) :: x
      if (self%scale <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
         x = quiet_nan()
      else if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else
         x = self%mu+self%scale*log(p/(1.0_dp-p))
      end if
   end function logistic_quantile

   function student_t_pdf(self, x) result(y)
      class(student_t_distribution), intent(in) :: self
      real(dp), intent(in) :: x
      real(dp) :: y, z, log_coef
      if (self%nu <= 0.0_dp .or. self%scale <= 0.0_dp) then
         y = quiet_nan()
      else
         z = (x-self%mu)/self%scale
         log_coef = log_gamma(0.5_dp*(self%nu+1.0_dp))- &
            log_gamma(0.5_dp*self%nu)-0.5_dp*log(self%nu*pi)-log(self%scale)
         y = exp(log_coef-0.5_dp*(self%nu+1.0_dp)*log(1.0_dp+z*z/self%nu))
      end if
   end function student_t_pdf

   function student_t_cdf(self, x) result(p)
      class(student_t_distribution), intent(in) :: self
      real(dp), intent(in) :: x
      real(dp) :: p
      if (self%nu <= 0.0_dp .or. self%scale <= 0.0_dp) then
         p = quiet_nan()
      else
         p = student_t_cdf_std((x-self%mu)/self%scale,self%nu)
      end if
   end function student_t_cdf

   function student_t_quantile(self, p) result(x)
      class(student_t_distribution), intent(in) :: self
      real(dp), intent(in) :: p
      real(dp) :: x
      real(dp) :: lo, hi, mid, pmid
      integer :: i

      if (self%nu <= 0.0_dp .or. self%scale <= 0.0_dp .or. &
          p < 0.0_dp .or. p > 1.0_dp) then
         x = quiet_nan()
         return
      else if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
         return
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
         return
      else if (abs(p-0.5_dp) <= epsilon(1.0_dp)) then
         x = self%mu
         return
      end if

      lo = -1.0_dp
      hi = 1.0_dp
      do while (student_t_cdf_std(lo,self%nu) > p)
         lo = 2.0_dp*lo
      end do
      do while (student_t_cdf_std(hi,self%nu) < p)
         hi = 2.0_dp*hi
      end do
      do i = 1, 160
         mid = 0.5_dp*(lo+hi)
         pmid = student_t_cdf_std(mid,self%nu)
         if (pmid < p) then
            lo = mid
         else
            hi = mid
         end if
         if (abs(hi-lo) <= 2.0e-13_dp*max(1.0_dp,abs(mid))) exit
      end do
      x = self%mu+self%scale*0.5_dp*(lo+hi)
   end function student_t_quantile

   subroutine initialize_callback_distribution(self, pdf_proc, cdf_proc, &
                                                quantile_proc, lower, upper)
      class(callback_distribution), intent(inout) :: self
      procedure(scalar_callback) :: pdf_proc, cdf_proc, quantile_proc
      real(dp), intent(in) :: lower, upper
      self%pdf_proc => pdf_proc
      self%cdf_proc => cdf_proc
      self%quantile_proc => quantile_proc
      self%lower = lower
      self%upper = upper
   end subroutine initialize_callback_distribution

   function callback_pdf(self, x) result(y)
      class(callback_distribution), intent(in) :: self
      real(dp), intent(in) :: x
      real(dp) :: y
      if (associated(self%pdf_proc)) then
         y = self%pdf_proc(x)
      else
         y = quiet_nan()
      end if
   end function callback_pdf

   function callback_cdf(self, x) result(p)
      class(callback_distribution), intent(in) :: self
      real(dp), intent(in) :: x
      real(dp) :: p
      if (associated(self%cdf_proc)) then
         p = self%cdf_proc(x)
      else
         p = quiet_nan()
      end if
   end function callback_cdf

   function callback_quantile(self, p) result(x)
      class(callback_distribution), intent(in) :: self
      real(dp), intent(in) :: p
      real(dp) :: x
      if (associated(self%quantile_proc)) then
         x = self%quantile_proc(p)
      else
         x = quiet_nan()
      end if
   end function callback_quantile

   function callback_lower(self) result(x)
      class(callback_distribution), intent(in) :: self
      real(dp) :: x
      x = self%lower
   end function callback_lower

   function callback_upper(self) result(x)
      class(callback_distribution), intent(in) :: self
      real(dp) :: x
      x = self%upper
   end function callback_upper

end module risk_distributions
