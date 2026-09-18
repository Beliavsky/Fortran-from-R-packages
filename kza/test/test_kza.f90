! SPDX-License-Identifier: GPL-3.0-only
program test_kza
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_quiet_nan, ieee_value
   use kza_api, only : dp, kz, kza, kzsv, kzft, kzs, kztp, periodogram, transfer_function, rlv
   implicit none

   real(dp), parameter :: tol = 1.0e-10_dp
   real(dp) :: x1(5)
   real(dp) :: xm(3)
   real(dp) :: a2(5, 5)
   real(dp) :: a3(5, 5, 1)
   real(dp) :: flat2(5, 5)
   real(dp), allocatable :: y1(:)
   real(dp), allocatable :: y2(:, :)
   real(dp), allocatable :: y3(:, :, :)
   real(dp), allocatable :: base2(:, :)
   real(dp), allocatable :: field(:, :)
   real(dp), allocatable :: baseline_field(:, :)
   real(dp), allocatable :: sym1(:, :)
   real(dp), allocatable :: sym2(:, :)
   real(dp), allocatable :: rotated(:, :)
   real(dp), allocatable :: sv(:)
   real(dp), allocatable :: smooth(:)
   real(dp), allocatable :: tf(:)
   real(dp), allocatable :: pg(:, :)
   real(dp), allocatable :: rv1(:)
   real(dp), allocatable :: rv2(:, :)
   real(dp), allocatable :: rv3(:, :, :)
   complex(dp), allocatable :: z(:)
   complex(dp), allocatable :: tp(:, :)
   integer :: i
   integer :: j

   x1 = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
   y1 = kz(x1, 3, k=1)
   call assert_close_vec(y1, [1.5_dp, 2.0_dp, 3.0_dp, 4.0_dp, 4.5_dp], tol, "kz 1D")

   xm = [1.0_dp, ieee_value(0.0_dp, ieee_quiet_nan), 3.0_dp]
   y1 = kz(xm, 3, k=1)
   call assert_close_vec(y1, [1.0_dp, 2.0_dp, 3.0_dp], tol, "kz missing-value averaging")

   do j = 1, 5
      do i = 1, 5
         a2(i, j) = real(i + 10 * j, dp)
      end do
   end do
   y2 = kz(a2, 3, k=1)
   call assert_close(y2(3, 3), sum(a2(2:4, 2:4)) / 9.0_dp, tol, "kz 2D center")

   a3(:, :, 1) = a2
   y3 = kz(a3, 3, k=1)
   call assert_close(maxval(abs(y3(:, :, 1) - y2)), 0.0_dp, tol, "kz single-slice 3D")

   y1 = kza(spread(3.5_dp, 1, 40), 7, k=3, min_size=1, impute_tails=.true.)
   call assert_close(maxval(abs(y1 - 3.5_dp)), 0.0_dp, tol, "kza constant 1D")

   flat2 = 0.0_dp
   base2 = kz(a2, 5, k=1)
   y2 = kza(a2, 5, y=flat2, k=1, min_size=0)
   call assert_close(maxval(abs(y2 - base2)), 0.0_dp, tol, "kza flat-baseline full window")

   y3 = kza(a3, 5, y=reshape(flat2, [5, 5, 1]), k=1, min_size=0)
   call assert_close(maxval(abs(y3(:, :, 1) - y2)), 0.0_dp, tol, "kza single-slice 3D")

   allocate(field(20, 20), baseline_field(20, 20))
   field = 0.0_dp
   field(6:15, 8:14) = 50.0_dp
   baseline_field = kz(field, 7, k=2)
   sym1 = kza(field, 7, y=baseline_field, k=2, min_size=1, symmetrize=.true.)
   rotated = rotate_test(field, 1)
   sym2 = kza(rotated, 7, y=rotate_test(baseline_field, 1), k=2, min_size=1, symmetrize=.true.)
   sym2 = rotate_test(sym2, 3)
   call assert_close(maxval(abs(sym1 - sym2)), 0.0_dp, 1.0e-9_dp, "kza rotational symmetrization")

   y2 = kza(field, 7, y=baseline_field, k=1, min_size=1)
   sym1 = kza(field, 7, y=baseline_field, k=1, min_size=1, normalize="max")
   call assert_close(maxval(abs(y2 - sym1)), 0.0_dp, tol, "kza default max normalization")

   baseline_field = 0.0_dp
   baseline_field(:, 11:20) = 1.0_dp
   baseline_field(4, 4) = 1000.0_dp
   y2 = kza(field, 7, y=baseline_field, k=1, min_size=1, normalize="max")
   sym1 = kza(field, 7, y=baseline_field, k=1, min_size=1, normalize="quantile")
   call assert_true(maxval(abs(y2 - sym1)) > 1.0e-12_dp, "kza quantile normalization changes dominated scale")

   deallocate(field, baseline_field)
   allocate(field(30, 30), baseline_field(30, 30))
   field = 0.0_dp
   field(14:17, 14:17) = 100.0_dp
   baseline_field = 0.0_dp
   baseline_field(15, 15) = 100.0_dp
   y2 = kza(field, 7, y=baseline_field, k=1, min_size=1, normalize="max")
   sym1 = kza(field, 7, y=baseline_field, k=1, min_size=1, normalize="quantile")
   call assert_close(maxval(abs(y2 - sym1)), 0.0_dp, tol, "kza sparse quantile fallback")

   y1 = kza(spread(2.0_dp, 1, 20), 5, k=1, min_size=1)
   call assert_true(all(ieee_is_nan(y1(1:5))), "kza leading tail marking")
   call assert_true(all(ieee_is_nan(y1(15:20))), "kza trailing tail marking")
   call assert_true(all(ieee_is_finite(y1(6:14))), "kza interior remains finite")

   sv = kzsv(spread(2.0_dp, 1, 20), spread(2.0_dp, 1, 20), 5, 1)
   call assert_close(maxval(abs(sv)), 0.0_dp, tol, "kzsv constant series")

   z = kzft(x1, f=0.0_dp, m=3.0_dp, k=1)
   call assert_close_vec(real(z, dp), [1.5_dp, 2.0_dp, 3.0_dp, 4.0_dp, 4.5_dp], tol, "kzft zero frequency")
   call assert_close(maxval(abs(aimag(z))), 0.0_dp, tol, "kzft zero-frequency imaginary part")

   smooth = kzs(x1, m=3.0_dp, k=1)
   call assert_close_vec(smooth, real(z, dp), tol, "kzs equals zero-frequency kzft")

   tf = transfer_function(5, 3, lamda=[0.25_dp], omega=0.25_dp)
   call assert_close(tf(1), 1.0_dp, tol, "transfer function at center frequency")

   pg = periodogram(spread(1.0_dp, 1, 8))
   call assert_close(pg(1, 1), 0.125_dp, tol, "periodogram first x-axis point")
   call assert_close(pg(1, 2), 8.0_dp, tol, "periodogram DC magnitude")
   call assert_close(maxval(abs(pg(2:, 2))), 0.0_dp, 1.0e-9_dp, "periodogram constant non-DC magnitude")

   tp = kztp([(sin(2.0_dp * acos(-1.0_dp) * real(i, dp) / 8.0_dp), i=1,24)], 4, 1)
   call assert_true(size(tp, 1) == 2 .and. size(tp, 2) == 2, "kztp default box shape")
   call assert_true(all(ieee_is_finite(real(tp, dp))), "kztp finite real values")
   call assert_true(all(ieee_is_finite(aimag(tp))), "kztp finite imaginary values")

   rv1 = rlv(spread(1.0_dp, 1, 7), 3)
   call assert_close(maxval(abs(rv1)), 0.0_dp, tol, "rlv clamp constant 1D")
   rv1 = rlv(spread(1.0_dp, 1, 7), 3, pad="zero")
   call assert_close(rv1(1), 1.0_dp / 3.0_dp, tol, "rlv zero-pad edge variance")
   call assert_close(rv1(4), 0.0_dp, tol, "rlv zero-pad interior variance")
   rv1 = rlv(spread(1.0_dp, 1, 3), 1)
   call assert_close(rv1(1), 0.5_dp, tol, "rlv krnl=1 left edge")
   call assert_close(rv1(2), 0.0_dp, tol, "rlv krnl=1 interior")

   rv2 = rlv(spread(spread(1.0_dp, 1, 4), 2, 4), 3)
   call assert_close(maxval(abs(rv2)), 0.0_dp, tol, "rlv clamp constant 2D")
   rv3 = rlv(spread(spread(spread(1.0_dp, 1, 3), 2, 3), 3, 2), 3)
   call assert_close(maxval(abs(rv3)), 0.0_dp, tol, "rlv clamp constant 3D")

   print '(a)', "All kza tests passed."

contains

   pure function rotate_test(a, rotations) result(b)
      real(dp), intent(in) :: a(:, :) !! Matrix rotated for the symmetrization regression test.
      integer, intent(in) :: rotations !! Number of 90-degree counterclockwise rotations modulo four.
      real(dp), allocatable :: b(:, :)
      integer :: i
      integer :: j
      integer :: r

      r = modulo(rotations, 4)
      select case (r)
      case (0)
         allocate(b(size(a, 1), size(a, 2)))
         b = a
      case (1)
         allocate(b(size(a, 2), size(a, 1)))
         do j = 1, size(b, 2)
            do i = 1, size(b, 1)
               b(i, j) = a(j, size(a, 2) - i + 1)
            end do
         end do
      case (2)
         allocate(b(size(a, 1), size(a, 2)))
         do j = 1, size(b, 2)
            do i = 1, size(b, 1)
               b(i, j) = a(size(a, 1) - i + 1, size(a, 2) - j + 1)
            end do
         end do
      case (3)
         allocate(b(size(a, 2), size(a, 1)))
         do j = 1, size(b, 2)
            do i = 1, size(b, 1)
               b(i, j) = a(size(a, 1) - j + 1, i)
            end do
         end do
      end select
   end function rotate_test

   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual !! Scalar value produced by the translated routine.
      real(dp), intent(in) :: expected !! Deterministic reference value expected by the test.
      real(dp), intent(in) :: tolerance !! Maximum permitted absolute error.
      character(len=*), intent(in) :: label !! Human-readable name identifying the tested behavior.

      if (.not. ieee_is_finite(actual) .or. abs(actual - expected) > tolerance) then
         write (*, '(a,2(1x,es24.16))') "FAIL " // trim(label) // ":", actual, expected
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_close_vec(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual(:) !! Vector produced by the translated routine.
      real(dp), intent(in) :: expected(:) !! Deterministic reference vector expected by the test.
      real(dp), intent(in) :: tolerance !! Maximum permitted elementwise absolute error.
      character(len=*), intent(in) :: label !! Human-readable name identifying the tested behavior.

      if (size(actual) /= size(expected)) then
         write (*, '(a)') "FAIL " // trim(label) // ": size mismatch"
         error stop 1
      end if
      if (any(.not. ieee_is_finite(actual)) .or. maxval(abs(actual - expected)) > tolerance) then
         write (*, '(a,1x,es24.16)') "FAIL " // trim(label) // ": max error", maxval(abs(actual - expected))
         error stop 1
      end if
   end subroutine assert_close_vec

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition !! Boolean condition that must be true for the test to pass.
      character(len=*), intent(in) :: label !! Human-readable name identifying the tested behavior.

      if (.not. condition) then
         write (*, '(a)') "FAIL " // trim(label)
         error stop 1
      end if
   end subroutine assert_true

end program test_kza
