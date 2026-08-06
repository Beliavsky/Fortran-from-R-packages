! Modern Fortran translation of the computational code in ADGofTest 0.3.
! The upstream package declares its license as "GPL" without a version.
! See LICENSE and licenses/ in this distribution.
module adgoftest
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_positive_inf
   implicit none
   private

   integer, parameter, public :: dp = kind(1.0d0)

   integer, parameter, public :: ad_success = 0
   integer, parameter, public :: ad_empty_sample = 1
   integer, parameter, public :: ad_probability_out_of_range = 2
   integer, parameter, public :: ad_nonfinite_input = 3
   integer, parameter, public :: ad_invalid_sample_size = 4

   type, public :: ad_test_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      integer :: n = 0
      integer :: status = ad_success
      character(len=160) :: message = ''
   end type ad_test_result

   abstract interface
      pure function scalar_cdf(x) result(probability)
         import :: dp
         real(dp), intent(in) :: x
         real(dp) :: probability
      end function scalar_cdf
   end interface

   interface ad_test
      module procedure ad_test_uniform
      module procedure ad_test_with_cdf
   end interface ad_test

   public :: ad_test
   public :: ad_test_uniform
   public :: ad_test_with_cdf
   public :: ad_statistic
   public :: ad_distribution_cdf
   public :: ad_p_value

contains

   function ad_statistic(probabilities, status, message, clip_probabilities) result(statistic)
      real(dp), intent(in) :: probabilities(:)
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message
      logical, intent(in), optional :: clip_probabilities
      real(dp) :: statistic

      real(dp), allocatable :: x(:)
      real(dp) :: lower
      real(dp) :: upper
      real(dp) :: term
      real(dp) :: epsilon
      integer :: i
      integer :: n
      integer :: local_status
      character(len=160) :: local_message
      logical :: do_clip

      n = size(probabilities)
      statistic = 0.0_dp
      local_status = ad_success
      local_message = ''
      do_clip = .false.
      if (present(clip_probabilities)) do_clip = clip_probabilities

      if (n == 0) then
         local_status = ad_empty_sample
         local_message = 'The sample must contain at least one observation.'
         call assign_status(status, message, local_status, local_message)
         return
      end if

      if (any(.not. ieee_is_finite(probabilities))) then
         local_status = ad_nonfinite_input
         local_message = 'All probabilities must be finite.'
         call assign_status(status, message, local_status, local_message)
         return
      end if

      if (any(probabilities < 0.0_dp) .or. any(probabilities > 1.0_dp)) then
         local_status = ad_probability_out_of_range
         local_message = 'All probabilities must lie in the closed interval [0, 1].'
         call assign_status(status, message, local_status, local_message)
         return
      end if

      allocate(x(n))
      x = probabilities
      call sort_ascending(x)

      if (do_clip) then
         epsilon = max(spacing(1.0_dp), tiny(1.0_dp))
         x = max(epsilon, min(1.0_dp - epsilon, x))
      else if (x(1) <= 0.0_dp .or. x(n) >= 1.0_dp) then
         statistic = ieee_value(0.0_dp, ieee_positive_inf)
         call assign_status(status, message, local_status, local_message)
         return
      end if

      statistic = 0.0_dp
      do i = 1, n
         lower = x(i)
         upper = x(n + 1 - i)
         term = real(2 * i - 1, dp) * log(lower * (1.0_dp - upper))
         statistic = statistic + term
      end do
      statistic = -statistic / real(n, dp) - real(n, dp)

      call assign_status(status, message, local_status, local_message)
   end function ad_statistic

   pure function ad_distribution_cdf(statistic, n, clamp_probability) result(cdf_value)
      real(dp), intent(in) :: statistic
      integer, intent(in) :: n
      logical, intent(in), optional :: clamp_probability
      real(dp) :: cdf_value

      real(dp) :: asymptotic_cdf
      real(dp) :: v
      real(dp) :: z
      logical :: do_clamp

      do_clamp = .false.
      if (present(clamp_probability)) do_clamp = clamp_probability

      if (n <= 0 .or. statistic <= 0.0_dp .or. .not. ieee_is_finite(statistic)) then
         cdf_value = 0.0_dp
         if (statistic >= huge(1.0_dp) / 2.0_dp) cdf_value = 1.0_dp
         return
      end if

      if (statistic < 2.0_dp) then
         asymptotic_cdf = exp(-1.2337141_dp / statistic) / sqrt(statistic) * &
            (2.00012_dp + (0.247105_dp - (0.0649821_dp - (0.0347962_dp - &
            (0.011672_dp - 0.00168691_dp * statistic) * statistic) * statistic) * &
            statistic) * statistic)
      else
         asymptotic_cdf = exp(-exp(1.0776_dp - (2.30695_dp - (0.43424_dp - &
            (0.082433_dp - (0.008056_dp - 0.0003146_dp * statistic) * statistic) * &
            statistic) * statistic) * statistic))
      end if

      if (asymptotic_cdf > 0.8_dp) then
         cdf_value = asymptotic_cdf + (-130.2137_dp + (745.2337_dp - &
            (1705.091_dp - (1950.646_dp - (1116.360_dp - 255.7844_dp * &
            asymptotic_cdf) * asymptotic_cdf) * asymptotic_cdf) * asymptotic_cdf) * &
            asymptotic_cdf) / real(n, dp)
      else
         z = 0.01265_dp + 0.1757_dp / real(n, dp)
         if (asymptotic_cdf < z) then
            v = asymptotic_cdf / z
            v = sqrt(v) * (1.0_dp - v) * (49.0_dp * v - 102.0_dp)
            cdf_value = asymptotic_cdf + v * (0.0037_dp / real(n * n, dp) + &
               0.00078_dp / real(n, dp) + 0.00006_dp) / real(n, dp)
         else
            v = (asymptotic_cdf - z) / (0.8_dp - z)
            v = -0.00022633_dp + (6.54034_dp - (14.6538_dp - (14.458_dp - &
               (8.259_dp - 1.91864_dp * v) * v) * v) * v) * v
            cdf_value = asymptotic_cdf + v * (0.04213_dp + 0.01365_dp / &
               real(n, dp)) / real(n, dp)
         end if
      end if

      if (do_clamp) cdf_value = max(0.0_dp, min(1.0_dp, cdf_value))
   end function ad_distribution_cdf

   pure function ad_p_value(statistic, n, clamp_probability) result(p_value)
      real(dp), intent(in) :: statistic
      integer, intent(in) :: n
      logical, intent(in), optional :: clamp_probability
      real(dp) :: p_value

      logical :: do_clamp

      do_clamp = .false.
      if (present(clamp_probability)) do_clamp = clamp_probability
      p_value = 1.0_dp - ad_distribution_cdf(statistic, n, do_clamp)
      if (do_clamp) p_value = max(0.0_dp, min(1.0_dp, p_value))
   end function ad_p_value

   subroutine ad_test_uniform(x, result, clip_probabilities, clamp_p_value)
      real(dp), intent(in) :: x(:)
      type(ad_test_result), intent(out) :: result
      logical, intent(in), optional :: clip_probabilities
      logical, intent(in), optional :: clamp_p_value

      logical :: do_clip
      logical :: do_clamp

      do_clip = .false.
      do_clamp = .false.
      if (present(clip_probabilities)) do_clip = clip_probabilities
      if (present(clamp_p_value)) do_clamp = clamp_p_value

      result%n = size(x)
      result%statistic = ad_statistic(x, result%status, result%message, do_clip)
      if (result%status /= ad_success) then
         result%p_value = 0.0_dp
         return
      end if
      result%p_value = ad_p_value(result%statistic, result%n, do_clamp)
   end subroutine ad_test_uniform

   subroutine ad_test_with_cdf(x, cdf, result, clip_probabilities, clamp_p_value)
      real(dp), intent(in) :: x(:)
      procedure(scalar_cdf) :: cdf
      type(ad_test_result), intent(out) :: result
      logical, intent(in), optional :: clip_probabilities
      logical, intent(in), optional :: clamp_p_value

      real(dp), allocatable :: probabilities(:)
      integer :: i
      logical :: do_clip
      logical :: do_clamp

      do_clip = .false.
      do_clamp = .false.
      if (present(clip_probabilities)) do_clip = clip_probabilities
      if (present(clamp_p_value)) do_clamp = clamp_p_value

      allocate(probabilities(size(x)))
      do i = 1, size(x)
         probabilities(i) = cdf(x(i))
      end do
      call ad_test_uniform(probabilities, result, do_clip, do_clamp)
   end subroutine ad_test_with_cdf

   subroutine sort_ascending(x)
      real(dp), intent(inout) :: x(:)

      real(dp) :: value
      integer :: i
      integer :: j

      do i = 2, size(x)
         value = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= value) exit
            x(j + 1) = x(j)
            j = j - 1
         end do
         x(j + 1) = value
      end do
   end subroutine sort_ascending

   subroutine assign_status(status, message, local_status, local_message)
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message
      integer, intent(in) :: local_status
      character(len=*), intent(in) :: local_message

      if (present(status)) status = local_status
      if (present(message)) message = local_message
   end subroutine assign_status

end module adgoftest
