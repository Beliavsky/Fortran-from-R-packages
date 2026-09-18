! SPDX-License-Identifier: GPL-3.0-only
module kza_rlv
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use kza_kinds, only : dp
   use kza_utils, only : quiet_nan, sample_variance_values
   implicit none
   private

   public :: rlv

   interface rlv
      module procedure rlv_1d
      module procedure rlv_2d
      module procedure rlv_3d
   end interface rlv

contains

   function rlv_1d(inpt, krnl, pad) result(variance)
      real(dp), intent(in) :: inpt(:) !! One-dimensional data whose rolling local sample variance is estimated.
      integer, intent(in) :: krnl !! Odd centered window length; krnl=1 selects the legacy one-sided two-cell mode.
      character(len=*), intent(in), optional :: pad !! Boundary rule, "clamp" by default or "zero" for legacy zero padding.
      real(dp), allocatable :: variance(:)
      integer :: i
      integer :: j
      integer :: ncount
      integer :: r
      real(dp) :: s1
      real(dp) :: s2
      logical :: zero_pad

      call validate_rlv_arguments(inpt, krnl, pad, zero_pad)
      allocate(variance(size(inpt)))
      if (krnl == 1) then
         do i = 1, size(inpt)
            variance(i) = legacy_one_sided_1d(inpt, i)
         end do
         return
      end if

      r = (krnl - 1) / 2
      do i = 1, size(inpt)
         s1 = 0.0_dp
         s2 = 0.0_dp
         do j = max(1, i - r), min(size(inpt), i + r)
            s1 = s1 + inpt(j)
            s2 = s2 + inpt(j) * inpt(j)
         end do
         if (zero_pad) then
            ncount = krnl
         else
            ncount = min(size(inpt), i + r) - max(1, i - r) + 1
         end if
         variance(i) = variance_from_sums(s1, s2, ncount)
      end do
   end function rlv_1d

   function rlv_2d(inpt, krnl, pad) result(variance)
      real(dp), intent(in) :: inpt(:, :) !! Two-dimensional data whose rolling local sample variance atlas is estimated.
      integer, intent(in) :: krnl !! Odd centered window length per dimension; krnl=1 selects the legacy four-corner mode.
      character(len=*), intent(in), optional :: pad !! Boundary rule, "clamp" by default or "zero" for legacy zero padding.
      real(dp), allocatable :: variance(:, :)
      integer :: i
      integer :: ii
      integer :: j
      integer :: jj
      integer :: ncount
      integer :: r
      real(dp) :: s1
      real(dp) :: s2
      logical :: zero_pad

      call validate_rlv_arguments(reshape(inpt, [size(inpt)]), krnl, pad, zero_pad)
      allocate(variance(size(inpt, 1), size(inpt, 2)))
      if (krnl == 1) then
         do j = 1, size(inpt, 2)
            do i = 1, size(inpt, 1)
               variance(i, j) = legacy_one_sided_2d(inpt, i, j)
            end do
         end do
         return
      end if

      r = (krnl - 1) / 2
      do j = 1, size(inpt, 2)
         do i = 1, size(inpt, 1)
            s1 = 0.0_dp
            s2 = 0.0_dp
            do jj = max(1, j - r), min(size(inpt, 2), j + r)
               do ii = max(1, i - r), min(size(inpt, 1), i + r)
                  s1 = s1 + inpt(ii, jj)
                  s2 = s2 + inpt(ii, jj) * inpt(ii, jj)
               end do
            end do
            if (zero_pad) then
               ncount = krnl ** 2
            else
               ncount = (min(size(inpt, 1), i + r) - max(1, i - r) + 1) * &
                  (min(size(inpt, 2), j + r) - max(1, j - r) + 1)
            end if
            variance(i, j) = variance_from_sums(s1, s2, ncount)
         end do
      end do
   end function rlv_2d

   function rlv_3d(inpt, krnl, pad) result(variance)
      real(dp), intent(in) :: inpt(:, :, :) !! Three-dimensional data whose rolling local sample variance atlas is estimated.
      integer, intent(in) :: krnl !! Odd centered window length per dimension; krnl=1 selects the legacy eight-corner mode.
      character(len=*), intent(in), optional :: pad !! Boundary rule, "clamp" by default or "zero" for legacy zero padding.
      real(dp), allocatable :: variance(:, :, :)
      integer :: i
      integer :: ii
      integer :: j
      integer :: jj
      integer :: l
      integer :: ll
      integer :: ncount
      integer :: r
      real(dp) :: s1
      real(dp) :: s2
      logical :: zero_pad

      call validate_rlv_arguments(reshape(inpt, [size(inpt)]), krnl, pad, zero_pad)
      allocate(variance(size(inpt, 1), size(inpt, 2), size(inpt, 3)))
      if (krnl == 1) then
         do l = 1, size(inpt, 3)
            do j = 1, size(inpt, 2)
               do i = 1, size(inpt, 1)
                  variance(i, j, l) = legacy_one_sided_3d(inpt, i, j, l)
               end do
            end do
         end do
         return
      end if

      r = (krnl - 1) / 2
      do l = 1, size(inpt, 3)
         do j = 1, size(inpt, 2)
            do i = 1, size(inpt, 1)
               s1 = 0.0_dp
               s2 = 0.0_dp
               do ll = max(1, l - r), min(size(inpt, 3), l + r)
                  do jj = max(1, j - r), min(size(inpt, 2), j + r)
                     do ii = max(1, i - r), min(size(inpt, 1), i + r)
                        s1 = s1 + inpt(ii, jj, ll)
                        s2 = s2 + inpt(ii, jj, ll) * inpt(ii, jj, ll)
                     end do
                  end do
               end do
               if (zero_pad) then
                  ncount = krnl ** 3
               else
                  ncount = (min(size(inpt, 1), i + r) - max(1, i - r) + 1) * &
                     (min(size(inpt, 2), j + r) - max(1, j - r) + 1) * &
                     (min(size(inpt, 3), l + r) - max(1, l - r) + 1)
               end if
               variance(i, j, l) = variance_from_sums(s1, s2, ncount)
            end do
         end do
      end do
   end function rlv_3d

   subroutine validate_rlv_arguments(flat, krnl, pad, zero_pad)
      real(dp), intent(in) :: flat(:) !! Flattened input used only to enforce upstream's clamp-mode missing-data restriction.
      integer, intent(in) :: krnl !! Requested rolling-window length.
      character(len=*), intent(in), optional :: pad !! Requested boundary policy, "clamp" or "zero".
      logical, intent(out) :: zero_pad !! True when the selected boundary policy is legacy zero padding.
      integer :: i
      character(len=:), allocatable :: mode

      if (krnl < 1) error stop "rlv: krnl must be at least 1"
      zero_pad = .false.
      if (present(pad)) then
         mode = trim(adjustl(pad))
         if (mode == "zero") then
            zero_pad = .true.
         else if (mode == "clamp") then
            zero_pad = .false.
         else
            error stop "rlv: pad must be 'clamp' or 'zero'"
         end if
      end if
      if (krnl == 1) return
      if (mod(krnl, 2) == 0) error stop "rlv: krnl must be odd"
      if (.not. zero_pad) then
         do i = 1, size(flat)
            if (ieee_is_nan(flat(i))) then
               error stop "rlv: pad='clamp' requires complete data"
            end if
         end do
      end if
   end subroutine validate_rlv_arguments

   pure real(dp) function variance_from_sums(s1, s2, ncount) result(v)
      real(dp), intent(in) :: s1 !! Sum of values in the effective window, excluding conceptual zero padding.
      real(dp), intent(in) :: s2 !! Sum of squared values in the effective window, excluding conceptual zero padding.
      integer, intent(in) :: ncount !! Sample size used in the variance denominator, including zero padding when requested.
      real(dp) :: raw

      if (ncount <= 1) then
         v = quiet_nan()
         return
      end if
      raw = (s2 - s1 * s1 / real(ncount, dp)) / real(ncount - 1, dp)
      if (ieee_is_nan(raw)) then
         v = raw
      else
         v = max(raw, 0.0_dp)
      end if
   end function variance_from_sums

   pure real(dp) function legacy_one_sided_1d(a, i) result(v)
      real(dp), intent(in) :: a(:) !! One-dimensional input accessed with conceptual zero padding at both ends.
      integer, intent(in) :: i !! Center position whose left and right two-cell sample variances are compared.
      real(dp) :: left_values(2)
      real(dp) :: right_values(2)
      real(dp) :: vleft
      real(dp) :: vright

      left_values = [padded_value_1d(a, i - 1), padded_value_1d(a, i)]
      right_values = [padded_value_1d(a, i), padded_value_1d(a, i + 1)]
      vleft = sample_variance_values(left_values)
      vright = sample_variance_values(right_values)
      if (ieee_is_nan(vleft) .or. ieee_is_nan(vright)) then
         v = quiet_nan()
      else
         v = max(vleft, vright)
      end if
   end function legacy_one_sided_1d

   pure real(dp) function legacy_one_sided_2d(a, i, j) result(v)
      real(dp), intent(in) :: a(:, :) !! Two-dimensional input accessed with conceptual zero padding around every edge.
      integer, intent(in) :: i !! First-dimension center index of the current cell.
      integer, intent(in) :: j !! Second-dimension center index of the current cell.
      real(dp) :: vals(4)
      real(dp) :: vv(4)
      integer :: di_low(4)
      integer :: dj_low(4)
      integer :: c
      integer :: di
      integer :: dj
      integer :: pos

      di_low = [-1, -1, 0, 0]
      dj_low = [0, -1, -1, 0]
      do c = 1, 4
         pos = 0
         do dj = dj_low(c), dj_low(c) + 1
            do di = di_low(c), di_low(c) + 1
               pos = pos + 1
               vals(pos) = padded_value_2d(a, i + di, j + dj)
            end do
         end do
         vv(c) = sample_variance_values(vals)
      end do
      if (any(ieee_is_nan(vv))) then
         v = quiet_nan()
      else
         v = maxval(vv)
      end if
   end function legacy_one_sided_2d

   pure real(dp) function legacy_one_sided_3d(a, i, j, l) result(v)
      real(dp), intent(in) :: a(:, :, :) !! Three-dimensional input accessed with conceptual zero padding around every edge.
      integer, intent(in) :: i !! First-dimension center index of the current voxel.
      integer, intent(in) :: j !! Second-dimension center index of the current voxel.
      integer, intent(in) :: l !! Third-dimension center index of the current voxel.
      real(dp) :: vals(8)
      real(dp) :: vv(8)
      integer :: c
      integer :: di
      integer :: dj
      integer :: dl
      integer :: low_i
      integer :: low_j
      integer :: low_l
      integer :: pos

      c = 0
      do low_l = -1, 0
         do low_j = -1, 0
            do low_i = -1, 0
               c = c + 1
               pos = 0
               do dl = low_l, low_l + 1
                  do dj = low_j, low_j + 1
                     do di = low_i, low_i + 1
                        pos = pos + 1
                        vals(pos) = padded_value_3d(a, i + di, j + dj, l + dl)
                     end do
                  end do
               end do
               vv(c) = sample_variance_values(vals)
            end do
         end do
      end do
      if (any(ieee_is_nan(vv))) then
         v = quiet_nan()
      else
         v = maxval(vv)
      end if
   end function legacy_one_sided_3d

   pure real(dp) function padded_value_1d(a, i) result(x)
      real(dp), intent(in) :: a(:) !! One-dimensional source data for a conceptual zero-padded lookup.
      integer, intent(in) :: i !! Requested source index, which may lie outside the actual array bounds.

      if (i < 1 .or. i > size(a)) then
         x = 0.0_dp
      else
         x = a(i)
      end if
   end function padded_value_1d

   pure real(dp) function padded_value_2d(a, i, j) result(x)
      real(dp), intent(in) :: a(:, :) !! Two-dimensional source data for a conceptual zero-padded lookup.
      integer, intent(in) :: i !! Requested first-dimension index, which may lie outside the actual array bounds.
      integer, intent(in) :: j !! Requested second-dimension index, which may lie outside the actual array bounds.

      if (i < 1 .or. i > size(a, 1) .or. j < 1 .or. j > size(a, 2)) then
         x = 0.0_dp
      else
         x = a(i, j)
      end if
   end function padded_value_2d

   pure real(dp) function padded_value_3d(a, i, j, l) result(x)
      real(dp), intent(in) :: a(:, :, :) !! Three-dimensional source data for a conceptual zero-padded lookup.
      integer, intent(in) :: i !! Requested first-dimension index, which may lie outside the actual array bounds.
      integer, intent(in) :: j !! Requested second-dimension index, which may lie outside the actual array bounds.
      integer, intent(in) :: l !! Requested third-dimension index, which may lie outside the actual array bounds.

      if (i < 1 .or. i > size(a, 1) .or. j < 1 .or. j > size(a, 2) .or. l < 1 .or. l > size(a, 3)) then
         x = 0.0_dp
      else
         x = a(i, j, l)
      end if
   end function padded_value_3d

end module kza_rlv
