! SPDX-License-Identifier: MIT
module r_rolling
   use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
   use r_kinds, only : dp
   use r_status, only : r_invalid_input, r_no_data, r_ok
   implicit none
   private

   public :: r_roll_mean_right, r_roll_mean_valid

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

end module r_rolling
