module rpart_utils
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
   use rpart_kinds, only : dp, i8
   implicit none
   private
   public :: finite_dp, nan_dp, argsort_real, sort_real, unique_sorted_real
   public :: weighted_mean, weighted_sse, validate_control, shuffle_int
   public :: cumulative_linear_interp, r_round_even

contains

   elemental logical function finite_dp(x)
      real(dp), intent(in) :: x
      finite_dp = ieee_is_finite(x)
   end function finite_dp

   elemental real(dp) function nan_dp()
      nan_dp = ieee_value(0.0_dp, ieee_quiet_nan)
   end function nan_dp

   subroutine argsort_real(x, idx)
      real(dp), intent(in) :: x(:)
      integer, allocatable, intent(out) :: idx(:)
      integer, allocatable :: tmp(:)
      integer :: n, i
      n = size(x)
      allocate(idx(n), tmp(n))
      if (n == 0) return
      idx = [(i, i=1,n)]
      call merge_sort_idx(x, idx, tmp, 1, n)
   end subroutine argsort_real

   recursive subroutine merge_sort_idx(x, idx, tmp, lo, hi)
      real(dp), intent(in) :: x(:)
      integer, intent(inout) :: idx(:), tmp(:)
      integer, intent(in) :: lo, hi
      integer :: mid, i, j, k
      if (lo >= hi) return
      mid = (lo + hi) / 2
      call merge_sort_idx(x, idx, tmp, lo, mid)
      call merge_sort_idx(x, idx, tmp, mid+1, hi)
      i = lo; j = mid + 1; k = lo
      do while (i <= mid .and. j <= hi)
         if (x(idx(i)) <= x(idx(j))) then
            tmp(k) = idx(i); i = i + 1
         else
            tmp(k) = idx(j); j = j + 1
         end if
         k = k + 1
      end do
      do while (i <= mid)
         tmp(k) = idx(i); i = i + 1; k = k + 1
      end do
      do while (j <= hi)
         tmp(k) = idx(j); j = j + 1; k = k + 1
      end do
      idx(lo:hi) = tmp(lo:hi)
   end subroutine merge_sort_idx

   subroutine sort_real(x)
      real(dp), intent(inout) :: x(:)
      integer, allocatable :: idx(:)
      real(dp), allocatable :: tmp(:)
      call argsort_real(x, idx)
      allocate(tmp(size(x)))
      if (size(x) > 0) tmp = x(idx)
      x = tmp
   end subroutine sort_real

   subroutine unique_sorted_real(x, out, reltol)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: out(:)
      real(dp), intent(in), optional :: reltol
      real(dp), allocatable :: z(:), temp(:)
      real(dp) :: tol, scale
      integer :: i, k
      tol = 0.0_dp
      if (present(reltol)) tol = max(0.0_dp, reltol)
      allocate(z(size(x)))
      z = x
      call sort_real(z)
      allocate(temp(size(z)))
      k = 0
      do i = 1, size(z)
         if (k == 0) then
            k = 1; temp(k) = z(i)
         else
            scale = max(1.0_dp, abs(temp(k)), abs(z(i)))
            if (abs(z(i)-temp(k)) > tol*scale) then
               k = k + 1; temp(k) = z(i)
            end if
         end if
      end do
      allocate(out(k))
      if (k > 0) out = temp(1:k)
   end subroutine unique_sorted_real

   real(dp) function weighted_mean(y, wt) result(ans)
      real(dp), intent(in) :: y(:), wt(:)
      real(dp) :: sw
      sw = sum(wt)
      if (sw > 0.0_dp) then
         ans = dot_product(y, wt) / sw
      else
         ans = 0.0_dp
      end if
   end function weighted_mean

   real(dp) function weighted_sse(y, wt, mean) result(ans)
      real(dp), intent(in) :: y(:), wt(:), mean
      ans = sum(wt * (y - mean)**2)
   end function weighted_sse

   subroutine validate_control(control, stat)
      use rpart_types, only : rpart_control
      type(rpart_control), intent(inout) :: control
      integer, intent(out), optional :: stat
      integer :: s
      s = 0
      if (control%maxcompete < 0) control%maxcompete = 0
      if (control%maxsurrogate < 0) control%maxsurrogate = 0
      if (control%xval < 0) control%xval = 0
      if (control%maxdepth < 1 .or. control%maxdepth > 30) s = 1
      if (control%minsplit < 1 .or. control%minbucket < 1) s = 2
      if (control%cp < 0.0_dp) s = 3
      if (control%usesurrogate < 0 .or. control%usesurrogate > 2) control%usesurrogate = 2
      if (control%surrogatestyle < 0 .or. control%surrogatestyle > 1) control%surrogatestyle = 0
      if (present(stat)) stat = s
   end subroutine validate_control

   subroutine shuffle_int(a, seed)
      integer, intent(inout) :: a(:)
      integer(i8), intent(inout) :: seed
      integer :: i, j, t
      real(dp) :: u
      do i = size(a), 2, -1
         call lcg_uniform(seed, u)
         j = 1 + int(u * real(i,dp))
         if (j > i) j = i
         t = a(i); a(i) = a(j); a(j) = t
      end do
   end subroutine shuffle_int

   subroutine lcg_uniform(seed, u)
      integer(i8), intent(inout) :: seed
      real(dp), intent(out) :: u
      integer(i8), parameter :: a = 2862933555777941757_i8
      integer(i8), parameter :: c = 3037000493_i8
      seed = a * seed + c
      u = real(iand(seed, int(z'7FFFFFFFFFFFFFFF',i8)), dp) / real(huge(1_i8), dp)
      if (u >= 1.0_dp) u = 0.9999999999999999_dp
   end subroutine lcg_uniform

   real(dp) function cumulative_linear_interp(x, y, x0) result(ans)
      real(dp), intent(in) :: x(:), y(:), x0
      integer :: i, n
      real(dp) :: f
      n = size(x)
      if (n == 0) then
         ans = 0.0_dp; return
      end if
      if (x0 <= x(1)) then
         ans = y(1); return
      end if
      if (x0 >= x(n)) then
         ans = y(n); return
      end if
      do i = 1, n-1
         if (x0 <= x(i+1)) then
            if (x(i+1) > x(i)) then
               f = (x0-x(i))/(x(i+1)-x(i))
               ans = y(i) + f*(y(i+1)-y(i))
            else
               ans = y(i)
            end if
            return
         end if
      end do
      ans = y(n)
   end function cumulative_linear_interp

   integer function r_round_even(x) result(k)
      real(dp), intent(in) :: x
      real(dp) :: fl, frac
      integer :: i
      fl = floor(x)
      i = int(fl)
      frac = x - fl
      if (frac < 0.5_dp) then
         k = i
      else if (frac > 0.5_dp) then
         k = i + 1
      else if (mod(i,2) == 0) then
         k = i
      else
         k = i + 1
      end if
   end function r_round_even

end module rpart_utils
