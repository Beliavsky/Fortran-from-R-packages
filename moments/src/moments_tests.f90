! SPDX-License-Identifier: GPL-2.0-or-later
module moments_tests
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
   use moments_kinds, only : dp
   use moments_probability, only : normal_cdf, normal_survival, chi_square_2_survival
   use moments_statistics, only : skewness, kurtosis, geary, moment
   use moments_status, only : MOMENTS_SUCCESS, MOMENTS_INVALID_ARGUMENT, &
      MOMENTS_INSUFFICIENT_DATA, MOMENTS_DEGENERATE_DATA, MOMENTS_NONFINITE_DATA
   implicit none
   private

   integer, parameter, public :: ALTERNATIVE_TWO_SIDED = 0
   integer, parameter, public :: ALTERNATIVE_LESS = 1
   integer, parameter, public :: ALTERNATIVE_GREATER = 2

   type, public :: moments_test_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: transformed = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      integer :: n = 0
      integer :: alternative = ALTERNATIVE_TWO_SIDED
      integer :: status = MOMENTS_SUCCESS
      character(len=48) :: method = ''
      character(len=64) :: alternative_text = ''
   end type moments_test_result

   public :: agostino_test, anscombe_test, bonett_test, jarque_test

contains

   real(dp) function quiet_nan() result(x)
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end function quiet_nan

   subroutine compact_finite(x, clean, status)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: clean(:)
      integer, intent(out) :: status
      integer :: i, count

      count = 0
      do i = 1, size(x)
         if (ieee_is_finite(x(i))) count = count + 1
      end do
      allocate(clean(count))
      count = 0
      do i = 1, size(x)
         if (ieee_is_finite(x(i))) then
            count = count + 1
            clean(count) = x(i)
         end if
      end do
      if (size(clean) == 0) then
         status = MOMENTS_NONFINITE_DATA
      else
         status = MOMENTS_SUCCESS
      end if
   end subroutine compact_finite

   subroutine set_tail_probability(result, z, positive_text, negative_text)
      type(moments_test_result), intent(inout) :: result
      real(dp), intent(in) :: z
      character(len=*), intent(in) :: positive_text, negative_text

      select case (result%alternative)
      case (ALTERNATIVE_TWO_SIDED)
         result%p_value = 2.0_dp * min(normal_cdf(z), normal_survival(z))
         result%alternative_text = 'two-sided difference'
      case (ALTERNATIVE_LESS)
         result%p_value = normal_survival(z)
         result%alternative_text = positive_text
      case (ALTERNATIVE_GREATER)
         result%p_value = normal_cdf(z)
         result%alternative_text = negative_text
      case default
         result%status = MOMENTS_INVALID_ARGUMENT
         result%p_value = quiet_nan()
      end select
   end subroutine set_tail_probability

   function agostino_test(x, alternative) result(result)
      real(dp), intent(in) :: x(:)
      integer, intent(in), optional :: alternative
      type(moments_test_result) :: result
      real(dp), allocatable :: clean(:)
      real(dp) :: s3, y, b2, w, d, a, z
      integer :: status

      result%method = "D'Agostino skewness test"
      if (present(alternative)) result%alternative = alternative
      call compact_finite(x, clean, status)
      result%n = size(clean)
      if (status /= MOMENTS_SUCCESS .or. result%n < 8 .or. result%n > 46340) then
         result%status = MOMENTS_INSUFFICIENT_DATA
         result%statistic = quiet_nan()
         result%transformed = quiet_nan()
         result%p_value = quiet_nan()
         return
      end if
      s3 = skewness(clean)
      if (.not. ieee_is_finite(s3)) then
         result%status = MOMENTS_DEGENERATE_DATA
         result%statistic = quiet_nan()
         result%transformed = quiet_nan()
         result%p_value = quiet_nan()
         return
      end if
      y = s3 * sqrt(real((result%n + 1) * (result%n + 3), dp) / &
         real(6 * (result%n - 2), dp))
      b2 = 3.0_dp * (real(result%n, dp)**2 + 27.0_dp * real(result%n, dp) - 70.0_dp) * &
         real(result%n + 1, dp) * real(result%n + 3, dp) / &
         (real(result%n - 2, dp) * real(result%n + 5, dp) * &
         real(result%n + 7, dp) * real(result%n + 9, dp))
      w = sqrt(-1.0_dp + sqrt(2.0_dp * (b2 - 1.0_dp)))
      d = 1.0_dp / sqrt(log(w))
      a = sqrt(2.0_dp / (w * w - 1.0_dp))
      z = d * asinh(y / a)
      result%statistic = s3
      result%transformed = z
      call set_tail_probability(result, z, 'positive skewness', 'negative skewness')
   end function agostino_test

   function anscombe_test(x, alternative) result(result)
      real(dp), intent(in) :: x(:)
      integer, intent(in), optional :: alternative
      type(moments_test_result) :: result
      real(dp), allocatable :: clean(:)
      real(dp) :: b, eb2, vb2, m3, a, xx, z
      integer :: n, status

      result%method = 'Anscombe-Glynn kurtosis test'
      if (present(alternative)) result%alternative = alternative
      call compact_finite(x, clean, status)
      n = size(clean)
      result%n = n
      if (status /= MOMENTS_SUCCESS .or. n < 4) then
         result%status = MOMENTS_INSUFFICIENT_DATA
         result%statistic = quiet_nan()
         result%transformed = quiet_nan()
         result%p_value = quiet_nan()
         return
      end if
      b = kurtosis(clean)
      if (.not. ieee_is_finite(b)) then
         result%status = MOMENTS_DEGENERATE_DATA
         result%statistic = quiet_nan()
         result%transformed = quiet_nan()
         result%p_value = quiet_nan()
         return
      end if
      eb2 = 3.0_dp * real(n - 1, dp) / real(n + 1, dp)
      vb2 = 24.0_dp * real(n, dp) * real(n - 2, dp) * real(n - 3, dp) / &
         (real(n + 1, dp)**2 * real(n + 3, dp) * real(n + 5, dp))
      m3 = 6.0_dp * (real(n, dp)**2 - 5.0_dp * real(n, dp) + 2.0_dp) / &
         (real(n + 7, dp) * real(n + 9, dp)) * &
         sqrt(6.0_dp * real(n + 3, dp) * real(n + 5, dp) / &
         (real(n, dp) * real(n - 2, dp) * real(n - 3, dp)))
      a = 6.0_dp + (8.0_dp / m3) * (2.0_dp / m3 + sqrt(1.0_dp + 4.0_dp / (m3 * m3)))
      xx = (b - eb2) / sqrt(vb2)
      z = (1.0_dp - 2.0_dp / (9.0_dp * a) - &
         ((1.0_dp - 2.0_dp / a) / (1.0_dp + xx * sqrt(2.0_dp / (a - 4.0_dp))))**(1.0_dp / 3.0_dp)) / &
         sqrt(2.0_dp / (9.0_dp * a))
      result%statistic = b
      result%transformed = z
      call set_tail_probability(result, z, 'kurtosis greater than 3', 'kurtosis lower than 3')
   end function anscombe_test

   function bonett_test(x, alternative) result(result)
      real(dp), intent(in) :: x(:)
      integer, intent(in), optional :: alternative
      type(moments_test_result) :: result
      real(dp), allocatable :: clean(:)
      real(dp) :: rho, tau, omega, z
      integer :: status

      result%method = 'Bonett-Seier test for Geary kurtosis'
      if (present(alternative)) result%alternative = alternative
      call compact_finite(x, clean, status)
      result%n = size(clean)
      if (status /= MOMENTS_SUCCESS .or. result%n < 2) then
         result%status = MOMENTS_INSUFFICIENT_DATA
         result%statistic = quiet_nan()
         result%transformed = quiet_nan()
         result%p_value = quiet_nan()
         return
      end if
      rho = sqrt(moment(clean, 2, central=.true.))
      tau = moment(clean, 1, central=.true., absolute=.true.)
      if (rho <= 0.0_dp .or. tau <= 0.0_dp) then
         result%status = MOMENTS_DEGENERATE_DATA
         result%statistic = quiet_nan()
         result%transformed = quiet_nan()
         result%p_value = quiet_nan()
         return
      end if
      omega = 13.29_dp * (log(rho) - log(tau))
      z = sqrt(real(result%n + 2, dp)) * (omega - 3.0_dp) / 3.54_dp
      result%statistic = tau
      result%transformed = z
      call set_tail_probability(result, z, 'Geary kurtosis greater than normal', &
         'Geary kurtosis lower than normal')
   end function bonett_test

   function jarque_test(x) result(result)
      real(dp), intent(in) :: x(:)
      type(moments_test_result) :: result
      real(dp) :: s, k, jb
      integer :: i

      result%method = 'Jarque-Bera normality test'
      result%alternative = ALTERNATIVE_GREATER
      result%alternative_text = 'greater'
      result%n = size(x)
      if (size(x) < 2) then
         result%status = MOMENTS_INSUFFICIENT_DATA
         result%statistic = quiet_nan()
         result%transformed = quiet_nan()
         result%p_value = quiet_nan()
         return
      end if
      do i = 1, size(x)
         if (.not. ieee_is_finite(x(i))) then
            result%status = MOMENTS_NONFINITE_DATA
            result%statistic = quiet_nan()
            result%transformed = quiet_nan()
            result%p_value = quiet_nan()
            return
         end if
      end do
      s = skewness(x)
      k = kurtosis(x)
      if (.not. ieee_is_finite(s) .or. .not. ieee_is_finite(k)) then
         result%status = MOMENTS_DEGENERATE_DATA
         result%statistic = quiet_nan()
         result%transformed = quiet_nan()
         result%p_value = quiet_nan()
         return
      end if
      jb = real(size(x), dp) / 6.0_dp * (s * s + 0.25_dp * (k - 3.0_dp)**2)
      result%statistic = jb
      result%transformed = quiet_nan()
      result%p_value = chi_square_2_survival(jb)
   end function jarque_test

end module moments_tests
