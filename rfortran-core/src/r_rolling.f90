! SPDX-License-Identifier: MIT
module r_rolling
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
   use, intrinsic :: ieee_arithmetic, only : ieee_negative_inf, ieee_positive_inf, ieee_quiet_nan, ieee_value
   use r_kinds, only : dp
   use r_descriptive, only : r_correlation, r_covariance, r_variance
   use r_status, only : r_invalid_input, r_no_data, r_ok
   implicit none
   private

   public :: r_roll_mean_right, r_roll_mean_valid
   public :: r_roll_max_right, r_roll_max_valid
   public :: r_roll_min_right, r_roll_min_valid
   public :: r_roll_sum_right, r_roll_sum_valid
   public :: r_roll_correlation_right, r_roll_correlation_valid
   public :: r_roll_covariance_right, r_roll_covariance_valid
   public :: r_roll_sd_right, r_roll_sd_valid
   public :: r_roll_variance_right, r_roll_variance_valid

contains

   subroutine r_roll_mean_valid(x, window, values, status)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: window
      real(dp), allocatable, intent(out) :: values(:)
      integer, intent(out), optional :: status
      integer :: i, output_size

      if (window < 1) then
         allocate(values(0))
         if (present(status)) status = r_invalid_input
         return
      end if
      if (window > size(x)) then
         allocate(values(0))
         if (present(status)) status = r_no_data
         return
      end if

      output_size = size(x) - window + 1
      allocate(values(output_size))
      do i = 1, output_size
         values(i) = sum(x(i:i + window - 1))/real(window, dp)
      end do
      if (present(status)) status = r_ok
   end subroutine r_roll_mean_valid

   subroutine r_roll_mean_right(x, window, values, status)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: window
      real(dp), allocatable, intent(out) :: values(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: valid(:)
      integer :: local_status

      allocate(values(size(x)))
      values = ieee_value(0.0_dp, ieee_quiet_nan)
      if (window < 1) then
         if (present(status)) status = r_invalid_input
         return
      end if
      if (window > size(x)) then
         if (present(status)) status = r_ok
         return
      end if

      call r_roll_mean_valid(x, window, valid, local_status)
      values(window:) = valid
      if (present(status)) status = local_status
   end subroutine r_roll_mean_right

   subroutine r_roll_sum_valid(x, window, values, na_rm, finite_only, status)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: window
      real(dp), allocatable, intent(out) :: values(:)
      logical, intent(in), optional :: na_rm, finite_only
      integer, intent(out), optional :: status
      integer :: i, output_size
      logical :: remove_na, require_finite

      remove_na = .false.
      if (present(na_rm)) remove_na = na_rm
      require_finite = .false.
      if (present(finite_only)) require_finite = finite_only
      if (window < 1) then
         allocate(values(0))
         if (present(status)) status = r_invalid_input
         return
      end if
      if (window > size(x)) then
         allocate(values(0))
         if (present(status)) status = r_no_data
         return
      end if

      output_size = size(x) - window + 1
      allocate(values(output_size))
      if (all(ieee_is_finite(x))) then
         values(1) = sum(x(:window))
         do i = 2, output_size
            values(i) = values(i - 1) - x(i - 1) + x(i + window - 1)
         end do
      else
         do i = 1, output_size
            values(i) = window_sum(x(i:i + window - 1), remove_na, require_finite)
         end do
      end if
      if (present(status)) status = r_ok
   end subroutine r_roll_sum_valid

   subroutine r_roll_sum_right(x, window, values, na_rm, finite_only, status)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: window
      real(dp), allocatable, intent(out) :: values(:)
      logical, intent(in), optional :: na_rm, finite_only
      integer, intent(out), optional :: status
      real(dp), allocatable :: valid(:)
      integer :: local_status

      allocate(values(size(x)))
      values = ieee_value(0.0_dp, ieee_quiet_nan)
      if (window < 1) then
         if (present(status)) status = r_invalid_input
         return
      end if
      if (window > size(x)) then
         if (present(status)) status = r_ok
         return
      end if
      call r_roll_sum_valid(x, window, valid, na_rm, finite_only, local_status)
      if (local_status == r_ok) values(window:) = valid
      if (present(status)) status = local_status
   end subroutine r_roll_sum_right

   subroutine r_roll_min_valid(x, window, values, na_rm, finite_only, status)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: window
      real(dp), allocatable, intent(out) :: values(:)
      logical, intent(in), optional :: na_rm, finite_only
      integer, intent(out), optional :: status
      integer :: i, output_size
      logical :: remove_na, require_finite

      remove_na = .false.
      if (present(na_rm)) remove_na = na_rm
      require_finite = .false.
      if (present(finite_only)) require_finite = finite_only
      if (window < 1) then
         allocate(values(0))
         if (present(status)) status = r_invalid_input
         return
      end if
      if (window > size(x)) then
         allocate(values(0))
         if (present(status)) status = r_no_data
         return
      end if

      output_size = size(x) - window + 1
      allocate(values(output_size))
      do i = 1, output_size
         values(i) = window_min(x(i:i + window - 1), remove_na, require_finite)
      end do
      if (present(status)) status = r_ok
   end subroutine r_roll_min_valid

   subroutine r_roll_min_right(x, window, values, na_rm, finite_only, status)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: window
      real(dp), allocatable, intent(out) :: values(:)
      logical, intent(in), optional :: na_rm, finite_only
      integer, intent(out), optional :: status
      real(dp), allocatable :: valid(:)
      integer :: local_status

      allocate(values(size(x)))
      values = ieee_value(0.0_dp, ieee_quiet_nan)
      if (window < 1) then
         if (present(status)) status = r_invalid_input
         return
      end if
      if (window > size(x)) then
         if (present(status)) status = r_ok
         return
      end if
      call r_roll_min_valid(x, window, valid, na_rm, finite_only, local_status)
      if (local_status == r_ok) values(window:) = valid
      if (present(status)) status = local_status
   end subroutine r_roll_min_right

   subroutine r_roll_max_valid(x, window, values, na_rm, finite_only, status)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: window
      real(dp), allocatable, intent(out) :: values(:)
      logical, intent(in), optional :: na_rm, finite_only
      integer, intent(out), optional :: status
      integer :: i, output_size
      logical :: remove_na, require_finite

      remove_na = .false.
      if (present(na_rm)) remove_na = na_rm
      require_finite = .false.
      if (present(finite_only)) require_finite = finite_only
      if (window < 1) then
         allocate(values(0))
         if (present(status)) status = r_invalid_input
         return
      end if
      if (window > size(x)) then
         allocate(values(0))
         if (present(status)) status = r_no_data
         return
      end if

      output_size = size(x) - window + 1
      allocate(values(output_size))
      do i = 1, output_size
         values(i) = window_max(x(i:i + window - 1), remove_na, require_finite)
      end do
      if (present(status)) status = r_ok
   end subroutine r_roll_max_valid

   subroutine r_roll_max_right(x, window, values, na_rm, finite_only, status)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: window
      real(dp), allocatable, intent(out) :: values(:)
      logical, intent(in), optional :: na_rm, finite_only
      integer, intent(out), optional :: status
      real(dp), allocatable :: valid(:)
      integer :: local_status

      allocate(values(size(x)))
      values = ieee_value(0.0_dp, ieee_quiet_nan)
      if (window < 1) then
         if (present(status)) status = r_invalid_input
         return
      end if
      if (window > size(x)) then
         if (present(status)) status = r_ok
         return
      end if
      call r_roll_max_valid(x, window, valid, na_rm, finite_only, local_status)
      if (local_status == r_ok) values(window:) = valid
      if (present(status)) status = local_status
   end subroutine r_roll_max_right

   subroutine r_roll_variance_valid(x, window, values, ddof, na_rm, status)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: window
      real(dp), allocatable, intent(out) :: values(:)
      integer, intent(in), optional :: ddof
      logical, intent(in), optional :: na_rm
      integer, intent(out), optional :: status
      integer :: degrees, i, output_size
      logical :: remove_na

      degrees = 1
      if (present(ddof)) degrees = ddof
      remove_na = .false.
      if (present(na_rm)) remove_na = na_rm
      if (window < 1 .or. degrees < 0) then
         allocate(values(0))
         if (present(status)) status = r_invalid_input
         return
      end if
      if (window > size(x)) then
         allocate(values(0))
         if (present(status)) status = r_no_data
         return
      end if

      output_size = size(x) - window + 1
      allocate(values(output_size))
      do i = 1, output_size
         values(i) = r_variance(x(i:i + window - 1), na_rm=remove_na, ddof=degrees)
      end do
      if (present(status)) status = r_ok
   end subroutine r_roll_variance_valid

   subroutine r_roll_variance_right(x, window, values, ddof, na_rm, status)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: window
      real(dp), allocatable, intent(out) :: values(:)
      integer, intent(in), optional :: ddof
      logical, intent(in), optional :: na_rm
      integer, intent(out), optional :: status
      real(dp), allocatable :: valid(:)
      integer :: local_status

      allocate(values(size(x)))
      values = ieee_value(0.0_dp, ieee_quiet_nan)
      if (window < 1) then
         if (present(status)) status = r_invalid_input
         return
      end if
      if (window > size(x)) then
         if (present(status)) status = r_ok
         return
      end if
      call r_roll_variance_valid(x, window, valid, ddof, na_rm, local_status)
      if (local_status == r_ok) values(window:) = valid
      if (present(status)) status = local_status
   end subroutine r_roll_variance_right

   subroutine r_roll_sd_valid(x, window, values, ddof, na_rm, status)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: window
      real(dp), allocatable, intent(out) :: values(:)
      integer, intent(in), optional :: ddof
      logical, intent(in), optional :: na_rm
      integer, intent(out), optional :: status
      integer :: local_status

      call r_roll_variance_valid(x, window, values, ddof, na_rm, local_status)
      where (values >= 0.0_dp) values = sqrt(values)
      if (present(status)) status = local_status
   end subroutine r_roll_sd_valid

   subroutine r_roll_sd_right(x, window, values, ddof, na_rm, status)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: window
      real(dp), allocatable, intent(out) :: values(:)
      integer, intent(in), optional :: ddof
      logical, intent(in), optional :: na_rm
      integer, intent(out), optional :: status
      integer :: local_status

      call r_roll_variance_right(x, window, values, ddof, na_rm, local_status)
      where (values >= 0.0_dp) values = sqrt(values)
      if (present(status)) status = local_status
   end subroutine r_roll_sd_right

   subroutine r_roll_covariance_valid(x, y, window, values, ddof, na_rm, finite_only, status)
      real(dp), intent(in) :: x(:), y(:)
      integer, intent(in) :: window
      real(dp), allocatable, intent(out) :: values(:)
      integer, intent(in), optional :: ddof
      logical, intent(in), optional :: na_rm, finite_only
      integer, intent(out), optional :: status
      integer :: degrees, i, output_size
      logical :: remove_na, require_finite

      degrees = 1
      if (present(ddof)) degrees = ddof
      remove_na = .false.
      if (present(na_rm)) remove_na = na_rm
      require_finite = .false.
      if (present(finite_only)) require_finite = finite_only
      if (size(x) /= size(y) .or. window < 1 .or. degrees < 0) then
         allocate(values(0))
         if (present(status)) status = r_invalid_input
         return
      end if
      if (window > size(x)) then
         allocate(values(0))
         if (present(status)) status = r_no_data
         return
      end if

      output_size = size(x) - window + 1
      allocate(values(output_size))
      do i = 1, output_size
         values(i) = r_covariance(x(i:i + window - 1), y(i:i + window - 1), remove_na, require_finite, degrees)
      end do
      if (present(status)) status = r_ok
   end subroutine r_roll_covariance_valid

   subroutine r_roll_covariance_right(x, y, window, values, ddof, na_rm, finite_only, status)
      real(dp), intent(in) :: x(:), y(:)
      integer, intent(in) :: window
      real(dp), allocatable, intent(out) :: values(:)
      integer, intent(in), optional :: ddof
      logical, intent(in), optional :: na_rm, finite_only
      integer, intent(out), optional :: status
      real(dp), allocatable :: valid(:)
      integer :: local_status

      if (size(x) /= size(y)) then
         allocate(values(0))
         if (present(status)) status = r_invalid_input
         return
      end if
      allocate(values(size(x)))
      values = ieee_value(0.0_dp, ieee_quiet_nan)
      if (window < 1) then
         if (present(status)) status = r_invalid_input
         return
      end if
      if (window > size(x)) then
         if (present(status)) status = r_ok
         return
      end if
      call r_roll_covariance_valid(x, y, window, valid, ddof, na_rm, finite_only, local_status)
      if (local_status == r_ok) values(window:) = valid
      if (present(status)) status = local_status
   end subroutine r_roll_covariance_right

   subroutine r_roll_correlation_valid(x, y, window, values, na_rm, finite_only, status)
      real(dp), intent(in) :: x(:), y(:)
      integer, intent(in) :: window
      real(dp), allocatable, intent(out) :: values(:)
      logical, intent(in), optional :: na_rm, finite_only
      integer, intent(out), optional :: status
      integer :: i, output_size
      logical :: remove_na, require_finite

      remove_na = .false.
      if (present(na_rm)) remove_na = na_rm
      require_finite = .false.
      if (present(finite_only)) require_finite = finite_only
      if (size(x) /= size(y) .or. window < 1) then
         allocate(values(0))
         if (present(status)) status = r_invalid_input
         return
      end if
      if (window > size(x)) then
         allocate(values(0))
         if (present(status)) status = r_no_data
         return
      end if

      output_size = size(x) - window + 1
      allocate(values(output_size))
      do i = 1, output_size
         values(i) = r_correlation(x(i:i + window - 1), y(i:i + window - 1), remove_na, require_finite)
      end do
      if (present(status)) status = r_ok
   end subroutine r_roll_correlation_valid

   subroutine r_roll_correlation_right(x, y, window, values, na_rm, finite_only, status)
      real(dp), intent(in) :: x(:), y(:)
      integer, intent(in) :: window
      real(dp), allocatable, intent(out) :: values(:)
      logical, intent(in), optional :: na_rm, finite_only
      integer, intent(out), optional :: status
      real(dp), allocatable :: valid(:)
      integer :: local_status

      if (size(x) /= size(y)) then
         allocate(values(0))
         if (present(status)) status = r_invalid_input
         return
      end if
      allocate(values(size(x)))
      values = ieee_value(0.0_dp, ieee_quiet_nan)
      if (window < 1) then
         if (present(status)) status = r_invalid_input
         return
      end if
      if (window > size(x)) then
         if (present(status)) status = r_ok
         return
      end if
      call r_roll_correlation_valid(x, y, window, valid, na_rm, finite_only, local_status)
      if (local_status == r_ok) values(window:) = valid
      if (present(status)) status = local_status
   end subroutine r_roll_correlation_right

   pure real(dp) function window_sum(x, remove_na, require_finite) result(value)
      real(dp), intent(in) :: x(:)
      logical, intent(in) :: remove_na, require_finite

      if (require_finite) then
         value = sum(x, mask=ieee_is_finite(x))
      else if (remove_na) then
         value = sum(x, mask=.not. ieee_is_nan(x))
      else if (any(ieee_is_nan(x))) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else
         value = sum(x)
      end if
   end function window_sum

   pure real(dp) function window_min(x, remove_na, require_finite) result(value)
      real(dp), intent(in) :: x(:)
      logical, intent(in) :: remove_na, require_finite
      logical :: keep(size(x))

      if (require_finite) then
         keep = ieee_is_finite(x)
      else if (remove_na) then
         keep = .not. ieee_is_nan(x)
      else if (any(ieee_is_nan(x))) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      else
         value = minval(x)
         return
      end if
      if (any(keep)) then
         value = minval(x, mask=keep)
      else
         value = ieee_value(0.0_dp, ieee_positive_inf)
      end if
   end function window_min

   pure real(dp) function window_max(x, remove_na, require_finite) result(value)
      real(dp), intent(in) :: x(:)
      logical, intent(in) :: remove_na, require_finite
      logical :: keep(size(x))

      if (require_finite) then
         keep = ieee_is_finite(x)
      else if (remove_na) then
         keep = .not. ieee_is_nan(x)
      else if (any(ieee_is_nan(x))) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      else
         value = maxval(x)
         return
      end if
      if (any(keep)) then
         value = maxval(x, mask=keep)
      else
         value = ieee_value(0.0_dp, ieee_negative_inf)
      end if
   end function window_max

end module r_rolling
