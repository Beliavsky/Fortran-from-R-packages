module pmwr_utils
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use pmwr_kinds, only : dp
   implicit none
   private
   public :: mean_value, sample_sd, copy_forward_rows, trunc_decimals
   public :: cumulative_product, finite_value, lower_ascii

contains

   pure real(dp) function mean_value(x) result(ans)
      real(dp), intent(in) :: x(:)
      if (size(x) == 0) then
         ans = 0.0_dp
      else
         ans = sum(x) / real(size(x), dp)
      end if
   end function mean_value

   pure real(dp) function sample_sd(x) result(ans)
      real(dp), intent(in) :: x(:)
      real(dp) :: m
      if (size(x) < 2) then
         ans = 0.0_dp
      else
         m = mean_value(x)
         ans = sqrt(sum((x - m)**2) / real(size(x) - 1, dp))
      end if
   end function sample_sd

   subroutine copy_forward_rows(x, valid)
      real(dp), intent(inout) :: x(:,:)
      logical, intent(in) :: valid(:)
      integer :: i
      if (size(valid) /= size(x, 1)) error stop "copy_forward_rows: size mismatch"
      do i = 2, size(x, 1)
         if (.not. valid(i)) x(i, :) = x(i - 1, :)
      end do
   end subroutine copy_forward_rows

   pure real(dp) function trunc_decimals(x, digits) result(ans)
      real(dp), intent(in) :: x
      integer, intent(in) :: digits
      real(dp) :: f
      f = 10.0_dp**real(digits, dp)
      ans = anint(aint(x * f) / f * f) / f
   end function trunc_decimals

   pure subroutine cumulative_product(x, ans)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: ans(:)
      integer :: i
      if (size(ans) /= size(x)) error stop "cumulative_product: size mismatch"
      if (size(x) == 0) return
      ans(1) = x(1)
      do i = 2, size(x)
         ans(i) = ans(i - 1) * x(i)
      end do
   end subroutine cumulative_product

   pure logical function finite_value(x) result(ok)
      real(dp), intent(in) :: x
      ok = ieee_is_finite(x)
   end function finite_value

   pure character(len=len(s)) function lower_ascii(s) result(ans)
      character(len=*), intent(in) :: s
      integer :: i, k
      ans = s
      do i = 1, len(s)
         k = iachar(ans(i:i))
         if (k >= iachar('A') .and. k <= iachar('Z')) ans(i:i) = achar(k + 32)
      end do
   end function lower_ascii

end module pmwr_utils
