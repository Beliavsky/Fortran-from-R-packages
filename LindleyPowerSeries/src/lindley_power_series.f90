! SPDX-License-Identifier: GPL-2.0-or-later
!
! Modern Fortran translation of the computational routines from
! LindleyPowerSeries 1.0.1 (GPL >= 2), by Saralees Nadarajah,
! Yuancheng Si, and Peihao Wang.
module lindley_power_series
   use lps_kinds, only : dp
   use lps_special, only : lambert_wm1, expm1_dp, log1p_dp, &
      log_expm1_pos, log_one_plus_p_expm1
   implicit none
   private

   public :: dp
   public :: lindley_cdf, lindley_pdf, lindley_quantile
   public :: plindleygeometric, dlindleygeometric, hlindleygeometric
   public :: qlindleygeometric, rlindleygeometric
   public :: plindleylogarithmic, dlindleylogarithmic, hlindleylogarithmic
   public :: qlindleylogarithmic, rlindleylogarithmic
   public :: plindleynb, dlindleynb, hlindleynb, qlindleynb, rlindleynb
   public :: plindleybinomial, dlindleybinomial, hlindleybinomial
   public :: qlindleybinomial, rlindleybinomial
   public :: plindleypoisson, dlindleypoisson, hlindleypoisson
   public :: qlindleypoisson, rlindleypoisson

contains

   pure elemental function lindley_cdf(x, lambda) result(p)
      real(dp), intent(in) :: x, lambda
      real(dp) :: p, log_surv

      if (x <= 0.0_dp) then
         p = 0.0_dp
         return
      end if
      log_surv = log1p_dp(lambda * x / (lambda + 1.0_dp)) - lambda * x
      p = -expm1_dp(log_surv)
      p = min(1.0_dp, max(0.0_dp, p))
   end function lindley_cdf

   pure elemental function lindley_pdf(x, lambda) result(f)
      real(dp), intent(in) :: x, lambda
      real(dp) :: f

      if (x < 0.0_dp) then
         f = 0.0_dp
      else
         f = lambda * lambda * (1.0_dp + x) * exp(-lambda*x) / &
             (lambda + 1.0_dp)
      end if
   end function lindley_pdf

   pure elemental function lindley_survival(x, lambda) result(s)
      real(dp), intent(in) :: x, lambda
      real(dp) :: s

      if (x <= 0.0_dp) then
         s = 1.0_dp
      else
         s = (1.0_dp + lambda*x/(lambda+1.0_dp)) * exp(-lambda*x)
      end if
   end function lindley_survival

   pure elemental function lindley_hazard(x, lambda) result(h)
      real(dp), intent(in) :: x, lambda
      real(dp) :: h

      if (x < 0.0_dp) then
         h = 0.0_dp
      else
         h = lambda*lambda*(1.0_dp+x) / (lambda+1.0_dp+lambda*x)
      end if
   end function lindley_hazard

   pure elemental function lindley_quantile(p, lambda) result(x)
      real(dp), intent(in) :: p, lambda
      real(dp) :: x, z, w

      if (p <= 0.0_dp) then
         x = 0.0_dp
         return
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
         return
      end if

      z = -(lambda + 1.0_dp) * exp(-(lambda + 1.0_dp)) * (1.0_dp - p)
      w = lambert_wm1(z)
      x = -w / lambda - 1.0_dp / lambda - 1.0_dp
      if (x < 0.0_dp .and. x > -64.0_dp*epsilon(1.0_dp)) x = 0.0_dp
   end function lindley_quantile

   pure elemental function plindleygeometric(x, lambda, theta, log_p) result(p)
      real(dp), intent(in) :: x, lambda, theta
      logical, intent(in), optional :: log_p
      real(dp) :: p, g
      logical :: lp

      g = lindley_cdf(x, lambda)
      p = g * (1.0_dp - theta) / (1.0_dp - theta*g)
      lp = .false.
      if (present(log_p)) lp = log_p
      if (lp) p = log(p)
   end function plindleygeometric

   pure elemental function dlindleygeometric(x, lambda, theta) result(f)
      real(dp), intent(in) :: x, lambda, theta
      real(dp) :: f, g, phi

      if (x < 0.0_dp) then
         f = 0.0_dp
         return
      end if
      g = lindley_cdf(x, lambda)
      phi = theta * g
      f = lindley_pdf(x, lambda) * (1.0_dp - theta) / (1.0_dp - phi)**2
   end function dlindleygeometric

   pure elemental function hlindleygeometric(x, lambda, theta) result(h)
      real(dp), intent(in) :: x, lambda, theta
      real(dp) :: h, g

      if (x < 0.0_dp) then
         h = 0.0_dp
         return
      end if
      g = lindley_cdf(x, lambda)
      h = lindley_hazard(x, lambda) * (1.0_dp - theta) / &
          (1.0_dp - theta*g)
   end function hlindleygeometric

   pure elemental function qlindleygeometric(p, lambda, theta) result(x)
      real(dp), intent(in) :: p, lambda, theta
      real(dp) :: x, u

      if (p <= 0.0_dp) then
         x = 0.0_dp
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else
         u = p / (1.0_dp - theta + p*theta)
         x = lindley_quantile(u, lambda)
      end if
   end function qlindleygeometric

   function rlindleygeometric(n, lambda, theta) result(x)
      integer, intent(in) :: n
      real(dp), intent(in) :: lambda, theta
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: u(:)

      allocate(x(n), u(n))
      call random_number(u)
      x = qlindleygeometric(u, lambda, theta)
   end function rlindleygeometric

   pure elemental function plindleylogarithmic(x, lambda, theta, log_p) result(p)
      real(dp), intent(in) :: x, lambda, theta
      logical, intent(in), optional :: log_p
      real(dp) :: p, g, aphi, atheta
      logical :: lp

      g = lindley_cdf(x, lambda)
      aphi = -log1p_dp(-theta*g)
      atheta = -log1p_dp(-theta)
      p = aphi / atheta
      lp = .false.
      if (present(log_p)) lp = log_p
      if (lp) p = log(p)
   end function plindleylogarithmic

   pure elemental function dlindleylogarithmic(x, lambda, theta) result(f)
      real(dp), intent(in) :: x, lambda, theta
      real(dp) :: f, g, atheta

      if (x < 0.0_dp) then
         f = 0.0_dp
         return
      end if
      g = lindley_cdf(x, lambda)
      atheta = -log1p_dp(-theta)
      f = theta * lindley_pdf(x, lambda) / &
          (atheta * (1.0_dp - theta*g))
   end function dlindleylogarithmic

   pure elemental function hlindleylogarithmic(x, lambda, theta) result(h)
      real(dp), intent(in) :: x, lambda, theta
      real(dp) :: h, g, s0, den, delta

      if (x < 0.0_dp) then
         h = 0.0_dp
         return
      end if
      g = lindley_cdf(x, lambda)
      s0 = lindley_survival(x, lambda)
      if (s0 <= tiny(1.0_dp)) then
         h = lindley_hazard(x, lambda)
         return
      end if
      den = 1.0_dp - theta*g
      delta = log1p_dp(theta*s0/(1.0_dp-theta))
      h = theta * lindley_hazard(x, lambda) * s0 / (den * delta)
   end function hlindleylogarithmic

   pure elemental function qlindleylogarithmic(p, lambda, theta) result(x)
      real(dp), intent(in) :: p, lambda, theta
      real(dp) :: x, u

      if (p <= 0.0_dp) then
         x = 0.0_dp
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else
         u = -expm1_dp(p * log1p_dp(-theta)) / theta
         x = lindley_quantile(u, lambda)
      end if
   end function qlindleylogarithmic

   function rlindleylogarithmic(n, lambda, theta) result(x)
      integer, intent(in) :: n
      real(dp), intent(in) :: lambda, theta
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: u(:)

      allocate(x(n), u(n))
      call random_number(u)
      x = qlindleylogarithmic(u, lambda, theta)
   end function rlindleylogarithmic

   pure elemental function plindleynb(x, lambda, theta, m, log_p) result(p)
      real(dp), intent(in) :: x, lambda, theta
      integer, intent(in) :: m
      logical, intent(in), optional :: log_p
      real(dp) :: p, g, ratio
      logical :: lp

      g = lindley_cdf(x, lambda)
      ratio = g * (1.0_dp - theta) / (1.0_dp - theta*g)
      p = ratio**m
      lp = .false.
      if (present(log_p)) lp = log_p
      if (lp) p = log(p)
   end function plindleynb

   pure elemental function dlindleynb(x, lambda, theta, m) result(f)
      real(dp), intent(in) :: x, lambda, theta
      integer, intent(in) :: m
      real(dp) :: f, g, den

      if (x < 0.0_dp) then
         f = 0.0_dp
         return
      end if
      g = lindley_cdf(x, lambda)
      den = 1.0_dp - theta*g
      if (m == 1) then
         f = lindley_pdf(x, lambda) * (1.0_dp - theta) / den**2
      else
         f = real(m, dp) * lindley_pdf(x, lambda) * g**(m-1) * &
             (1.0_dp - theta)**m / den**(m+1)
      end if
   end function dlindleynb

   pure elemental function hlindleynb(x, lambda, theta, m) result(h)
      real(dp), intent(in) :: x, lambda, theta
      integer, intent(in) :: m
      real(dp) :: h, g, s0, den, log_r, cdf, surv

      if (x < 0.0_dp) then
         h = 0.0_dp
         return
      end if
      g = lindley_cdf(x, lambda)
      s0 = lindley_survival(x, lambda)
      if (s0 <= tiny(1.0_dp)) then
         h = lindley_hazard(x, lambda)
         return
      end if
      den = 1.0_dp - theta*g
      log_r = log1p_dp(-s0/den)
      cdf = exp(real(m,dp)*log_r)
      surv = -expm1_dp(real(m,dp)*log_r)
      if (g <= tiny(1.0_dp)) then
         h = dlindleynb(x,lambda,theta,m)
      else
         h = real(m,dp) * lindley_hazard(x,lambda) * s0 * cdf / &
             (g * den * surv)
      end if
   end function hlindleynb

   pure elemental function qlindleynb(p, lambda, theta, m) result(x)
      real(dp), intent(in) :: p, lambda, theta
      integer, intent(in) :: m
      real(dp) :: x, r, u

      if (p <= 0.0_dp) then
         x = 0.0_dp
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else
         r = exp(log(p) / real(m, dp))
         u = r / (1.0_dp - theta + r*theta)
         x = lindley_quantile(u, lambda)
      end if
   end function qlindleynb

   function rlindleynb(n, lambda, theta, m) result(x)
      integer, intent(in) :: n, m
      real(dp), intent(in) :: lambda, theta
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: u(:)

      allocate(x(n), u(n))
      call random_number(u)
      x = qlindleynb(u, lambda, theta, m)
   end function rlindleynb

   pure elemental function plindleybinomial(x, lambda, theta, m, log_p) result(p)
      real(dp), intent(in) :: x, lambda, theta
      integer, intent(in) :: m
      logical, intent(in), optional :: log_p
      real(dp) :: p, g, log_num, log_den
      logical :: lp

      g = lindley_cdf(x, lambda)
      if (g <= 0.0_dp) then
         p = 0.0_dp
      else
         log_num = log_expm1_pos(real(m,dp) * log1p_dp(theta*g))
         log_den = log_expm1_pos(real(m,dp) * log1p_dp(theta))
         p = exp(log_num - log_den)
      end if
      lp = .false.
      if (present(log_p)) lp = log_p
      if (lp) p = log(p)
   end function plindleybinomial

   pure elemental function dlindleybinomial(x, lambda, theta, m) result(f)
      real(dp), intent(in) :: x, lambda, theta
      integer, intent(in) :: m
      real(dp) :: f, g, log_atheta

      if (x < 0.0_dp) then
         f = 0.0_dp
         return
      end if
      g = lindley_cdf(x, lambda)
      log_atheta = log_expm1_pos(real(m,dp) * log1p_dp(theta))
      f = theta * real(m,dp) * lindley_pdf(x, lambda) * &
          exp(real(m-1,dp)*log1p_dp(theta*g) - log_atheta)
   end function dlindleybinomial

   pure elemental function hlindleybinomial(x, lambda, theta, m) result(h)
      real(dp), intent(in) :: x, lambda, theta
      integer, intent(in) :: m
      real(dp) :: h, g, s0, den, delta

      if (x < 0.0_dp) then
         h = 0.0_dp
         return
      end if
      g = lindley_cdf(x, lambda)
      s0 = lindley_survival(x, lambda)
      if (s0 <= tiny(1.0_dp)) then
         h = lindley_hazard(x, lambda)
         return
      end if
      den = 1.0_dp + theta*g
      delta = real(m,dp) * log1p_dp(theta*s0/den)
      h = theta * real(m,dp) * lindley_hazard(x,lambda) * s0 / &
          (den * expm1_dp(delta))
   end function hlindleybinomial

   pure elemental function qlindleybinomial(p, lambda, theta, m) result(x)
      real(dp), intent(in) :: p, lambda, theta
      integer, intent(in) :: m
      real(dp) :: x, ltheta, log_one_plus_y, u

      if (p <= 0.0_dp) then
         x = 0.0_dp
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else
         ltheta = real(m,dp) * log1p_dp(theta)
         log_one_plus_y = log_one_plus_p_expm1(p, ltheta)
         u = expm1_dp(log_one_plus_y / real(m,dp)) / theta
         x = lindley_quantile(u, lambda)
      end if
   end function qlindleybinomial

   function rlindleybinomial(n, lambda, theta, m) result(x)
      integer, intent(in) :: n, m
      real(dp), intent(in) :: lambda, theta
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: u(:)

      allocate(x(n), u(n))
      call random_number(u)
      x = qlindleybinomial(u, lambda, theta, m)
   end function rlindleybinomial

   pure elemental function plindleypoisson(x, lambda, theta, log_p) result(p)
      real(dp), intent(in) :: x, lambda, theta
      logical, intent(in), optional :: log_p
      real(dp) :: p, g, phi, log_num, log_den
      logical :: lp

      g = lindley_cdf(x, lambda)
      if (g <= 0.0_dp) then
         p = 0.0_dp
      else
         phi = theta * g
         log_num = log_expm1_pos(phi)
         log_den = log_expm1_pos(theta)
         p = exp(log_num - log_den)
      end if
      lp = .false.
      if (present(log_p)) lp = log_p
      if (lp) p = log(p)
   end function plindleypoisson

   pure elemental function dlindleypoisson(x, lambda, theta) result(f)
      real(dp), intent(in) :: x, lambda, theta
      real(dp) :: f, g, log_den

      if (x < 0.0_dp) then
         f = 0.0_dp
         return
      end if
      g = lindley_cdf(x, lambda)
      log_den = log_expm1_pos(theta)
      f = theta * lindley_pdf(x, lambda) * exp(theta*g - log_den)
   end function dlindleypoisson

   pure elemental function hlindleypoisson(x, lambda, theta) result(h)
      real(dp), intent(in) :: x, lambda, theta
      real(dp) :: h, s0, delta

      if (x < 0.0_dp) then
         h = 0.0_dp
         return
      end if
      s0 = lindley_survival(x, lambda)
      if (s0 <= tiny(1.0_dp)) then
         h = lindley_hazard(x, lambda)
         return
      end if
      delta = theta * s0
      h = theta * lindley_hazard(x,lambda) * s0 / expm1_dp(delta)
   end function hlindleypoisson

   pure elemental function qlindleypoisson(p, lambda, theta) result(x)
      real(dp), intent(in) :: p, lambda, theta
      real(dp) :: x, u

      if (p <= 0.0_dp) then
         x = 0.0_dp
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else
         u = log_one_plus_p_expm1(p, theta) / theta
         x = lindley_quantile(u, lambda)
      end if
   end function qlindleypoisson

   function rlindleypoisson(n, lambda, theta) result(x)
      integer, intent(in) :: n
      real(dp), intent(in) :: lambda, theta
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: u(:)

      allocate(x(n), u(n))
      call random_number(u)
      x = qlindleypoisson(u, lambda, theta)
   end function rlindleypoisson

end module lindley_power_series
