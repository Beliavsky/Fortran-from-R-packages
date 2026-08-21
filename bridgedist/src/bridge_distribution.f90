! SPDX-License-Identifier: GPL-2.0-or-later
!
! Modern Fortran translation of the computational routines in the R package
! bridgedist 0.1.3 by Bruce Swihart.
module bridge_distribution
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf, &
      ieee_is_finite, ieee_is_nan
   use bridgedist_kinds, only : dp
   implicit none
   private

   real(dp), parameter :: pi = acos(-1.0_dp)

   public :: dbridge, pbridge, qbridge, rbridge
   public :: dbridge_recycle, pbridge_recycle, qbridge_recycle
   public :: bridge_mean, bridge_variance

   interface rbridge
      module procedure rbridge_scalar_phi
      module procedure rbridge_vector_phi
   end interface rbridge

contains

   elemental real(dp) function dbridge(x, phi, log_value) result(ans)
      real(dp), intent(in) :: x
      real(dp), intent(in), optional :: phi
      logical, intent(in), optional :: log_value
      real(dp) :: ph, a, t, logden
      logical :: want_log

      ph = 0.5_dp
      if (present(phi)) ph = phi
      want_log = .false.
      if (present(log_value)) want_log = log_value

      if (.not. valid_phi(ph)) then
         ans = nan_value()
         return
      end if

      if (.not. ieee_is_finite(x)) then
         if (ieee_is_nan(x)) then
            ans = nan_value()
         else if (want_log) then
            ans = ieee_value(ans, ieee_negative_inf)
         else
            ans = 0.0_dp
         end if
         return
      end if

      a = pi * ph
      t = ph * x
      logden = log_cosh_plus_cos(t, cos(a))
      ans = log(sin(a)) - log(2.0_dp * pi) - logden
      if (.not. want_log) ans = exp(ans)
   end function dbridge

   elemental real(dp) function pbridge(q, phi, lower_tail, log_p) result(ans)
      real(dp), intent(in) :: q
      real(dp), intent(in), optional :: phi
      logical, intent(in), optional :: lower_tail, log_p
      real(dp) :: ph, prob
      logical :: lower, want_log

      ph = 0.5_dp
      if (present(phi)) ph = phi
      lower = .true.
      if (present(lower_tail)) lower = lower_tail
      want_log = .false.
      if (present(log_p)) want_log = log_p

      if (.not. valid_phi(ph)) then
         ans = nan_value()
         return
      end if

      if (ieee_is_nan(q)) then
         ans = nan_value()
         return
      end if

      if (lower) then
         prob = bridge_cdf_lower(q, ph)
      else
         ! The bridge distribution is symmetric about zero.  Evaluating the
         ! upper tail as F(-q) avoids cancellation from 1-F(q).
         prob = bridge_cdf_lower(-q, ph)
      end if

      if (want_log) then
         if (prob <= 0.0_dp) then
            ans = ieee_value(ans, ieee_negative_inf)
         else
            ans = log(prob)
         end if
      else
         ans = prob
      end if
   end function pbridge

   elemental real(dp) function qbridge(p, phi, lower_tail, log_p) result(ans)
      real(dp), intent(in) :: p
      real(dp), intent(in), optional :: phi
      logical, intent(in), optional :: lower_tail, log_p
      real(dp) :: ph, prob
      logical :: lower, input_log

      ph = 0.5_dp
      if (present(phi)) ph = phi
      lower = .true.
      if (present(lower_tail)) lower = lower_tail
      input_log = .false.
      if (present(log_p)) input_log = log_p

      if (.not. valid_phi(ph)) then
         ans = nan_value()
         return
      end if

      if (input_log) then
         if (p > 0.0_dp .or. ieee_is_nan(p)) then
            ans = nan_value()
            return
         end if
         prob = exp(p)
      else
         prob = p
      end if

      if (prob < 0.0_dp .or. prob > 1.0_dp .or. ieee_is_nan(prob)) then
         ans = nan_value()
         return
      end if

      ans = bridge_quantile_lower(prob, ph)
      if (.not. lower) ans = -ans
   end function qbridge

   subroutine rbridge_scalar_phi(x, phi)
      real(dp), intent(out) :: x(:)
      real(dp), intent(in), optional :: phi
      real(dp) :: ph
      real(dp), allocatable :: u(:)

      ph = 0.5_dp
      if (present(phi)) ph = phi
      if (.not. valid_phi(ph)) then
         x = nan_value()
         return
      end if

      allocate(u(size(x)))
      call random_number(u)
      x = qbridge(u, ph)
   end subroutine rbridge_scalar_phi

   subroutine rbridge_vector_phi(x, phi)
      real(dp), intent(out) :: x(:)
      real(dp), intent(in) :: phi(:)
      real(dp) :: u
      integer :: i, j

      if (size(phi) == 0) then
         x = nan_value()
         return
      end if

      do i = 1, size(x)
         j = 1 + modulo(i - 1, size(phi))
         call random_number(u)
         x(i) = qbridge(u, phi(j))
      end do
   end subroutine rbridge_vector_phi

   function dbridge_recycle(x, phi, log_value) result(ans)
      real(dp), intent(in) :: x(:), phi(:)
      logical, intent(in), optional :: log_value
      real(dp), allocatable :: ans(:)
      integer :: i, n
      logical :: lv

      lv = .false.
      if (present(log_value)) lv = log_value
      if (size(x) == 0 .or. size(phi) == 0) then
         allocate(ans(0))
         return
      end if
      n = max(size(x), size(phi))
      allocate(ans(n))
      do i = 1, n
         ans(i) = dbridge(x(1 + modulo(i - 1, size(x))), phi(1 + modulo(i - 1, size(phi))), lv)
      end do
   end function dbridge_recycle

   function pbridge_recycle(q, phi, lower_tail, log_p) result(ans)
      real(dp), intent(in) :: q(:), phi(:)
      logical, intent(in), optional :: lower_tail, log_p
      real(dp), allocatable :: ans(:)
      integer :: i, n
      logical :: lower, lp

      lower = .true.
      if (present(lower_tail)) lower = lower_tail
      lp = .false.
      if (present(log_p)) lp = log_p
      if (size(q) == 0 .or. size(phi) == 0) then
         allocate(ans(0))
         return
      end if
      n = max(size(q), size(phi))
      allocate(ans(n))
      do i = 1, n
         ans(i) = pbridge(q(1 + modulo(i - 1, size(q))), phi(1 + modulo(i - 1, size(phi))), lower, lp)
      end do
   end function pbridge_recycle

   function qbridge_recycle(p, phi, lower_tail, log_p) result(ans)
      real(dp), intent(in) :: p(:), phi(:)
      logical, intent(in), optional :: lower_tail, log_p
      real(dp), allocatable :: ans(:)
      integer :: i, n
      logical :: lower, lp

      lower = .true.
      if (present(lower_tail)) lower = lower_tail
      lp = .false.
      if (present(log_p)) lp = log_p
      if (size(p) == 0 .or. size(phi) == 0) then
         allocate(ans(0))
         return
      end if
      n = max(size(p), size(phi))
      allocate(ans(n))
      do i = 1, n
         ans(i) = qbridge(p(1 + modulo(i - 1, size(p))), phi(1 + modulo(i - 1, size(phi))), lower, lp)
      end do
   end function qbridge_recycle

   elemental real(dp) function bridge_mean(phi) result(ans)
      real(dp), intent(in), optional :: phi
      real(dp) :: ph

      ph = 0.5_dp
      if (present(phi)) ph = phi
      if (valid_phi(ph)) then
         ans = 0.0_dp
      else
         ans = nan_value()
      end if
   end function bridge_mean

   elemental real(dp) function bridge_variance(phi) result(ans)
      real(dp), intent(in), optional :: phi
      real(dp) :: ph

      ph = 0.5_dp
      if (present(phi)) ph = phi
      if (valid_phi(ph)) then
         ans = pi**2 * (1.0_dp / ph**2 - 1.0_dp) / 3.0_dp
      else
         ans = nan_value()
      end if
   end function bridge_variance

   elemental logical function valid_phi(phi) result(ok)
      real(dp), intent(in) :: phi
      ok = ieee_is_finite(phi) .and. phi > 0.0_dp .and. phi < 1.0_dp
   end function valid_phi

   elemental real(dp) function bridge_cdf_lower(x, phi) result(prob)
      real(dp), intent(in) :: x, phi
      real(dp) :: a, t, e, theta

      if (.not. ieee_is_finite(x)) then
         if (ieee_is_nan(x)) then
            prob = nan_value()
         else if (x < 0.0_dp) then
            prob = 0.0_dp
         else
            prob = 1.0_dp
         end if
         return
      end if

      a = pi * phi
      t = phi * x
      if (t >= 0.0_dp) then
         e = exp(-t)
         theta = atan2(sin(a), e + cos(a))
      else
         e = exp(t)
         theta = atan2(e * sin(a), 1.0_dp + e * cos(a))
      end if
      prob = theta / a
      prob = min(1.0_dp, max(0.0_dp, prob))
   end function bridge_cdf_lower

   elemental real(dp) function bridge_quantile_lower(prob, phi) result(q)
      real(dp), intent(in) :: prob, phi
      real(dp) :: a

      if (prob <= 0.0_dp) then
         q = ieee_value(q, ieee_negative_inf)
         return
      else if (prob >= 1.0_dp) then
         q = ieee_value(q, ieee_positive_inf)
         return
      end if

      a = pi * phi
      q = (log(sin(a * prob)) - log(sin(a * (1.0_dp - prob)))) / phi
   end function bridge_quantile_lower

   elemental real(dp) function log_cosh_plus_cos(t, c) result(ans)
      real(dp), intent(in) :: t, c
      real(dp) :: at, e

      at = abs(t)
      if (at < 20.0_dp) then
         ans = log(cosh(t) + c)
      else
         e = exp(-at)
         ans = at - log(2.0_dp) + log(1.0_dp + 2.0_dp * c * e + e * e)
      end if
   end function log_cosh_plus_cos

   elemental real(dp) function nan_value() result(x)
      x = ieee_value(x, ieee_quiet_nan)
   end function nan_value

end module bridge_distribution
