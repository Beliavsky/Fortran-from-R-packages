! SPDX-License-Identifier: MIT
module cgnm_utils
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use cgnm_kinds, only : dp
   implicit none
   private
   public :: seed_rng, finite_real, median_value, quantile_value, sort_index
   public :: column_mean_sd, elbow_index, uniform01
contains
   subroutine seed_rng(seed)
      integer, intent(in) :: seed
      integer :: n, i
      integer, allocatable :: put(:)
      call random_seed(size=n)
      allocate(put(n))
      do i = 1, n
         put(i) = modulo(seed + 104729*i + 37*i*i, huge(1)-1)
         if (put(i) <= 0) put(i) = i + 1
      end do
      call random_seed(put=put)
   end subroutine seed_rng

   real(dp) function uniform01() result(u)
      call random_number(u)
   end function uniform01

   logical function finite_real(x) result(ok)
      real(dp), intent(in) :: x
      ok = ieee_is_finite(x)
   end function finite_real

   subroutine sort_index(x, idx)
      real(dp), intent(in) :: x(:)
      integer, intent(out) :: idx(size(x))
      integer :: i, j, k, t
      do i = 1, size(x)
         idx(i) = i
      end do
      do i = 1, size(x)-1
         k = i
         do j = i+1, size(x)
            if (x(idx(j)) < x(idx(k))) k = j
         end do
         if (k /= i) then
            t = idx(i); idx(i) = idx(k); idx(k) = t
         end if
      end do
   end subroutine sort_index

   real(dp) function median_value(x) result(v)
      real(dp), intent(in) :: x(:)
      integer, allocatable :: idx(:)
      integer :: n
      n = size(x)
      if (n == 0) then
         v = 0.0_dp
         return
      end if
      allocate(idx(n)); call sort_index(x, idx)
      if (mod(n,2) == 1) then
         v = x(idx((n+1)/2))
      else
         v = 0.5_dp*(x(idx(n/2)) + x(idx(n/2+1)))
      end if
   end function median_value

   real(dp) function quantile_value(x, prob) result(v)
      real(dp), intent(in) :: x(:), prob
      integer, allocatable :: idx(:)
      real(dp) :: h, frac
      integer :: n, j
      n = size(x)
      if (n == 0) then
         v = 0.0_dp; return
      end if
      allocate(idx(n)); call sort_index(x, idx)
      if (prob <= 0.0_dp) then
         v = x(idx(1)); return
      else if (prob >= 1.0_dp) then
         v = x(idx(n)); return
      end if
      ! R quantile default type 7
      h = 1.0_dp + real(n-1,dp)*prob
      j = int(floor(h))
      frac = h-real(j,dp)
      if (j >= n) then
         v = x(idx(n))
      else
         v = (1.0_dp-frac)*x(idx(j)) + frac*x(idx(j+1))
      end if
   end function quantile_value

   subroutine column_mean_sd(x, meanv, sdv)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(out) :: meanv(size(x,2)), sdv(size(x,2))
      integer :: j, n
      real(dp) :: s
      n = size(x,1)
      do j = 1, size(x,2)
         meanv(j) = sum(x(:,j))/real(n,dp)
         if (n > 1) then
            s = sum((x(:,j)-meanv(j))**2)/real(n-1,dp)
            sdv(j) = sqrt(max(0.0_dp,s))
         else
            sdv(j) = 0.0_dp
         end if
      end do
   end subroutine column_mean_sd

   integer function elbow_index(v) result(ind)
      real(dp), intent(in) :: v(:)
      real(dp) :: best, area
      integer :: i, n
      n = size(v)
      if (n <= 1) then
         ind = 1; return
      end if
      best = huge(1.0_dp); ind = 1
      do i = 1, n-1
         area = (v(1)+v(i))*real(i,dp) + &
                (v(i+1)+v(n))*real(n-i,dp)
         if (area < best) then
            best = area; ind = i
         end if
      end do
   end function elbow_index
end module cgnm_utils
