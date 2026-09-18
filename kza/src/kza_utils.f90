! SPDX-License-Identifier: GPL-3.0-only
module kza_utils
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_quiet_nan, ieee_value
   use kza_kinds, only : dp
   implicit none
   private

   public :: adaptive_factor
   public :: complex_is_nan
   public :: mean_finite_1d
   public :: normalizer
   public :: quiet_nan
   public :: r_round_even
   public :: sample_variance_values

contains

   pure elemental logical function complex_is_nan(z) result(is_nan)
      complex(dp), intent(in) :: z !! Complex value tested using R-like NA/NaN semantics for KZFT.

      is_nan = ieee_is_nan(real(z, dp)) .or. ieee_is_nan(aimag(z))
   end function complex_is_nan

   pure elemental real(dp) function quiet_nan() result(x)
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end function quiet_nan

   pure elemental real(dp) function adaptive_factor(d, scale) result(a)
      real(dp), intent(in) :: d !! Local absolute difference metric used to shrink an adaptive window.
      real(dp), intent(in) :: scale !! Positive global difference normalizer; nonpositive values mean no detected structure.

      if (scale <= 0.0_dp) then
         a = 1.0_dp
      else
         a = 1.0_dp - d / scale
      end if
   end function adaptive_factor

   pure elemental integer function r_round_even(x) result(ir)
      real(dp), intent(in) :: x !! Real value rounded to the nearest integer with R's ties-to-even convention.
      real(dp) :: ax
      real(dp) :: frac
      integer :: base
      integer :: sgn

      if (x < 0.0_dp) then
         sgn = -1
         ax = -x
      else
         sgn = 1
         ax = x
      end if

      base = floor(ax)
      frac = ax - real(base, dp)
      if (frac < 0.5_dp) then
         ir = base
      else if (frac > 0.5_dp) then
         ir = base + 1
      else if (mod(base, 2) == 0) then
         ir = base
      else
         ir = base + 1
      end if
      ir = sgn * ir
   end function r_round_even

   pure real(dp) function mean_finite_1d(x) result(avg)
      real(dp), intent(in) :: x(:) !! Values to average; non-finite entries are skipped like the compiled KZ moving averages.
      integer :: i
      integer :: nfinite
      real(dp) :: total

      total = 0.0_dp
      nfinite = 0
      do i = 1, size(x)
         if (ieee_is_finite(x(i))) then
            total = total + x(i)
            nfinite = nfinite + 1
         end if
      end do
      if (nfinite == 0) then
         avg = quiet_nan()
      else
         avg = total / real(nfinite, dp)
      end if
   end function mean_finite_1d

   pure real(dp) function sample_variance_values(x) result(v)
      real(dp), intent(in) :: x(:) !! Values whose ordinary sample variance is required; NaNs propagate.
      integer :: i
      real(dp) :: avg
      real(dp) :: ss

      if (size(x) <= 1) then
         v = quiet_nan()
         return
      end if
      avg = sum(x) / real(size(x), dp)
      ss = 0.0_dp
      do i = 1, size(x)
         ss = ss + (x(i) - avg) ** 2
      end do
      v = ss / real(size(x) - 1, dp)
   end function sample_variance_values

   pure real(dp) function normalizer(d, prob) result(scale)
      real(dp), intent(in) :: d(:) !! Difference metrics; only finite entries participate in the normalizer.
      real(dp), intent(in) :: prob !! Quantile probability in (0,1), or a value outside that interval to request the maximum.
      real(dp), allocatable :: finite_values(:)
      real(dp) :: h
      real(dp) :: q
      real(dp) :: max_value
      integer :: i
      integer :: lo
      integer :: nf

      nf = count(ieee_is_finite(d))
      if (nf == 0) then
         scale = quiet_nan()
         return
      end if

      max_value = -huge(1.0_dp)
      do i = 1, size(d)
         if (ieee_is_finite(d(i))) max_value = max(max_value, d(i))
      end do

      if (.not. (prob > 0.0_dp .and. prob < 1.0_dp)) then
         scale = max_value
         return
      end if

      allocate(finite_values(nf))
      nf = 0
      do i = 1, size(d)
         if (ieee_is_finite(d(i))) then
            nf = nf + 1
            finite_values(nf) = d(i)
         end if
      end do
      call quicksort_real(finite_values, 1, size(finite_values))

      if (size(finite_values) == 1) then
         q = finite_values(1)
      else
         h = real(size(finite_values) - 1, dp) * prob
         lo = floor(h) + 1
         if (lo >= size(finite_values)) then
            q = finite_values(size(finite_values))
         else
            q = finite_values(lo) + (h - floor(h)) * (finite_values(lo + 1) - finite_values(lo))
         end if
      end if

      if (q <= 0.0_dp) then
         scale = max_value
      else
         scale = q
      end if
   end function normalizer

   pure recursive subroutine quicksort_real(a, left, right)
      real(dp), intent(inout) :: a(:) !! Real work array sorted in ascending order in place.
      integer, intent(in) :: left !! First array index in the current quicksort partition.
      integer, intent(in) :: right !! Last array index in the current quicksort partition.
      integer :: i
      integer :: j
      real(dp) :: pivot
      real(dp) :: tmp

      if (left >= right) return
      i = left
      j = right
      pivot = a((left + right) / 2)
      do
         do while (a(i) < pivot)
            i = i + 1
         end do
         do while (a(j) > pivot)
            j = j - 1
         end do
         if (i <= j) then
            tmp = a(i)
            a(i) = a(j)
            a(j) = tmp
            i = i + 1
            j = j - 1
         end if
         if (i > j) exit
      end do
      if (left < j) call quicksort_real(a, left, j)
      if (i < right) call quicksort_real(a, i, right)
   end subroutine quicksort_real

end module kza_utils
