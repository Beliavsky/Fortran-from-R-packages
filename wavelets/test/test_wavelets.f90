! SPDX-License-Identifier: GPL-2.0-or-later
program test_wavelets
   use wavelets
   implicit none

   integer :: failures

   failures = 0
   call test_filter_api(failures)
   call test_extension(failures)
   call test_low_level_dwt(failures)
   call test_low_level_modwt(failures)
   call test_dwt_roundtrip(failures)
   call test_modwt_roundtrip(failures)
   call test_alignment(failures)
   call test_mra(failures)

   if (failures /= 0) then
      write (*, '(a,i0)') 'wavelets tests failed: ', failures
      error stop 1
   end if
   write (*, '(a)') 'wavelets tests passed'

contains

   subroutine test_filter_api(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions across the test program.

      type(wt_filter_type) :: filter
      type(wt_filter_type) :: mod_filter
      type(wt_filter_type) :: equivalent
      real(dp), allocatable :: qmf(:)
      real(dp), allocatable :: shifts(:)
      character(len=4), parameter :: names(26) = [character(len=4) :: &
         "haar", "d2", "d4", "d6", "d8", "d10", "d12", "d14", "d16", "d18", "d20", &
         "la8", "la10", "la12", "la14", "la16", "la18", "la20", "bl14", "bl18", "bl20", &
         "c6", "c12", "c18", "c24", "c30"]
      integer, parameter :: lengths(26) = [2, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, &
                                           8, 10, 12, 14, 16, 18, 20, 14, 18, 20, 6, 12, 18, 24, 30]
      integer :: shift
      integer :: ierr
      integer :: k

      do k = 1, size(names)
         call wt_filter_named(trim(names(k)), filter, ierr=ierr)
         call assert_true(ierr == 0, "named filter construction: " // trim(names(k)), failures)
         call assert_true(filter%l == lengths(k), "named filter length: " // trim(names(k)), failures)
         call assert_true(abs(sum(filter%g**2) - 1.0_dp) < 1.0e-9_dp, &
                          "named filter scaling norm: " // trim(names(k)), failures)
         call assert_true(abs(sum(filter%h**2) - 1.0_dp) < 1.0e-9_dp, &
                          "named filter wavelet norm: " // trim(names(k)), failures)
      end do

      call wt_filter_named('haar', filter, ierr=ierr)
      call assert_true(ierr == 0, 'haar filter status', failures)
      call assert_true(filter%l == 2, 'haar length', failures)
      call assert_close_vector(filter%g, [0.7071067811865475_dp, 0.7071067811865475_dp], &
                               1.0e-15_dp, 'haar scaling coefficients', failures)
      call assert_close_vector(filter%h, [0.7071067811865475_dp, -0.7071067811865475_dp], &
                               1.0e-15_dp, 'haar wavelet coefficients', failures)

      call wt_filter_qmf([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], qmf)
      call assert_close_vector(qmf, [-4.0_dp, 3.0_dp, -2.0_dp, 1.0_dp], 0.0_dp, 'qmf convention', failures)

      call wt_filter_named('la8', filter, ierr=ierr)
      call assert_true(ierr == 0, 'la8 filter status', failures)
      call wt_filter_equivalent(filter, 3, equivalent, ierr)
      call assert_true(ierr == 0, 'equivalent filter status', failures)
      call assert_true(equivalent%l == 50, 'la8 level-3 equivalent length', failures)
      call assert_true(equivalent%level == 3, 'equivalent filter level metadata', failures)

      call wt_filter_shift(filter, [1, 2, 3], shifts, wavelet=.true., ierr=ierr)
      call assert_true(ierr == 0, 'dwt shift status', failures)
      call assert_close_vector(shifts, [2.0_dp, 2.0_dp, 3.0_dp], 0.0_dp, 'la8 dwt shifts', failures)

      call wt_filter_named('la8', mod_filter, modwt=.true., ierr=ierr)
      call wt_filter_shift(mod_filter, [1, 2, 3], shifts, wavelet=.true., modwt=.true., ierr=ierr)
      call assert_close_vector(shifts, [4.0_dp, 11.0_dp, 25.0_dp], 0.0_dp, 'la8 modwt shifts', failures)

      call waveletshift_dwt(8, 2, shift, ierr=ierr)
      call assert_true(ierr == 0 .and. shift == 11, 'waveletshift.dwt la8 level 2', failures)
      call scalingshift_dwt(8, 2, shift, ierr=ierr)
      call assert_true(ierr == 0 .and. shift == 9, 'scalingshift.dwt la8 level 2', failures)
   end subroutine test_filter_api

   subroutine test_extension(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions across the test program.

      real(dp) :: x(3, 1)
      real(dp), allocatable :: y(:, :)
      integer :: ierr

      x(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp]

      call extend_series(x, y, method='reflection', ierr=ierr)
      call assert_true(ierr == 0, 'reflection extension status', failures)
      call assert_close_vector(y(:, 1), [1.0_dp, 2.0_dp, 3.0_dp, 3.0_dp, 2.0_dp, 1.0_dp], &
                               0.0_dp, 'reflection extension', failures)

      call extend_series(x, y, method='reflection.inverse', ierr=ierr)
      call assert_close_vector(y(:, 1), [1.0_dp, 2.0_dp, 3.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], &
                               0.0_dp, 'inverse reflection extension', failures)

      call extend_series(x, y, method='periodic', length_mode='arbitrary', n=8, ierr=ierr)
      call assert_close_vector(y(:, 1), [1.0_dp, 2.0_dp, 3.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 1.0_dp, 2.0_dp], &
                               0.0_dp, 'periodic arbitrary extension', failures)

      call extend_series(x, y, method='zeros', length_mode='arbitrary', n=5, ierr=ierr)
      call assert_close_vector(y(:, 1), [1.0_dp, 2.0_dp, 3.0_dp, 0.0_dp, 0.0_dp], &
                               0.0_dp, 'zero extension', failures)

      call extend_series(x, y, method='mean', length_mode='arbitrary', n=5, ierr=ierr)
      call assert_close_vector(y(:, 1), [1.0_dp, 2.0_dp, 3.0_dp, 2.0_dp, 2.0_dp], &
                               0.0_dp, 'mean extension', failures)

      call extend_series(x, y, method='reflection', length_mode='powerof2', j=2, ierr=ierr)
      call assert_close_vector(y(:, 1), [1.0_dp, 2.0_dp, 3.0_dp, 3.0_dp], &
                               0.0_dp, 'power-of-two extension', failures)
   end subroutine test_extension

   subroutine test_low_level_dwt(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions across the test program.

      type(wt_filter_type) :: filter
      real(dp) :: x(8)
      real(dp), allocatable :: w(:)
      real(dp), allocatable :: v(:)
      real(dp), allocatable :: reconstructed(:)
      integer :: ierr
      integer :: i

      do i = 1, 8
         x(i) = real(i, dp)
      end do
      call wt_filter_named('haar', filter, ierr=ierr)
      call dwt_forward(x, filter, w, v, ierr)
      call assert_true(ierr == 0, 'dwt.forward status', failures)
      call assert_close_vector(w, [0.7071067811865475_dp, 0.7071067811865475_dp, &
                                   0.7071067811865475_dp, 0.7071067811865475_dp], &
                               2.0e-15_dp, 'haar dwt.forward wavelets', failures)
      call assert_close_vector(v, [2.1213203435596424_dp, 4.9497474683058327_dp, &
                                   7.7781745930520225_dp, 10.606601717798213_dp], &
                               2.0e-14_dp, 'haar dwt.forward scaling', failures)
      call dwt_backward(w, v, filter, reconstructed, ierr)
      call assert_true(ierr == 0, 'dwt.backward status', failures)
      call assert_close_vector(reconstructed, x, 2.0e-14_dp, 'haar low-level roundtrip', failures)

      call wt_filter_named('haar', filter, modwt=.true., ierr=ierr)
      call dwt_forward(x, filter, w, v, ierr)
      call assert_true(ierr == 2, 'dwt rejects modwt-normalized filter', failures)
   end subroutine test_low_level_dwt

   subroutine test_low_level_modwt(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions across the test program.

      type(wt_filter_type) :: filter
      real(dp) :: x(8)
      real(dp), allocatable :: w(:)
      real(dp), allocatable :: v(:)
      real(dp), allocatable :: reconstructed(:)
      integer :: ierr
      integer :: i

      do i = 1, 8
         x(i) = real(i, dp)
      end do
      call wt_filter_named('haar', filter, modwt=.true., ierr=ierr)
      call modwt_forward(x, filter, 2, w, v, ierr)
      call assert_true(ierr == 0, 'modwt.forward status', failures)
      call assert_close_vector(w, [-3.0_dp, -3.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
                               2.0e-15_dp, 'haar modwt.forward wavelets', failures)
      call assert_close_vector(v, [4.0_dp, 5.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp], &
                               2.0e-15_dp, 'haar modwt.forward scaling', failures)
      call modwt_backward(w, v, filter, 2, reconstructed, ierr)
      call assert_true(ierr == 0, 'modwt.backward status', failures)
      call assert_close_vector(reconstructed, x, 2.0e-14_dp, 'haar low-level modwt roundtrip', failures)

      call wt_filter_named('haar', filter, ierr=ierr)
      call modwt_forward(x, filter, 2, w, v, ierr)
      call assert_true(ierr == 2, 'modwt rejects dwt-normalized filter', failures)
   end subroutine test_low_level_modwt

   subroutine test_dwt_roundtrip(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions across the test program.

      type(wavelet_transform_type) :: wt
      real(dp) :: x(64, 2)
      real(dp), allocatable :: reconstructed(:, :)
      integer :: ierr
      integer :: i

      do i = 1, size(x, 1)
         x(i, 1) = real(i, dp) / 8.0_dp
         x(i, 2) = real(mod(3 * i, 17), dp) / 16.0_dp
      end do

      call dwt_named(x, wt, filter_name='la8', boundary='periodic', ierr=ierr)
      call assert_true(ierr == 0 .and. wt%level == 3, 'default dwt level selection', failures)
      call dwt_named(x, wt, filter_name='la8', n_levels=3, boundary='periodic', ierr=ierr)
      call assert_true(ierr == 0, 'periodic dwt status', failures)
      call idwt(wt, reconstructed, ierr)
      call assert_true(ierr == 0, 'periodic idwt status', failures)
      call assert_close_matrix(reconstructed, x, 1.0e-12_dp, 'periodic dwt/idwt roundtrip', failures)

      call dwt_named(x, wt, filter_name='la8', n_levels=3, boundary='reflection', ierr=ierr)
      call assert_true(ierr == 0, 'reflection dwt status', failures)
      call idwt(wt, reconstructed, ierr)
      call assert_true(ierr == 0, 'reflection idwt status', failures)
      call assert_close_matrix(reconstructed, x, 1.0e-12_dp, 'reflection dwt/idwt roundtrip', failures)
   end subroutine test_dwt_roundtrip

   subroutine test_modwt_roundtrip(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions across the test program.

      type(wavelet_transform_type) :: wt
      real(dp) :: x(32, 2)
      real(dp), allocatable :: reconstructed(:, :)
      integer :: ierr
      integer :: i

      do i = 1, size(x, 1)
         x(i, 1) = real(i, dp) / 10.0_dp
         x(i, 2) = real(mod(5 * i, 19), dp) / 20.0_dp
      end do

      call modwt_named(x, wt, filter_name='la8', boundary='periodic', ierr=ierr)
      call assert_true(ierr == 0 .and. wt%level == 2, 'default modwt level selection', failures)
      call modwt_named(x, wt, filter_name='la8', n_levels=3, boundary='periodic', ierr=ierr)
      call assert_true(ierr == 0, 'periodic modwt status', failures)
      call imodwt(wt, reconstructed, ierr)
      call assert_true(ierr == 0, 'periodic imodwt status', failures)
      call assert_close_matrix(reconstructed, x, 1.0e-12_dp, 'periodic modwt/imodwt roundtrip', failures)

      call modwt_named(x, wt, filter_name='la8', n_levels=3, boundary='reflection', ierr=ierr)
      call assert_true(ierr == 0, 'reflection modwt status', failures)
      call imodwt(wt, reconstructed, ierr)
      call assert_true(ierr == 0, 'reflection imodwt status', failures)
      call assert_close_matrix(reconstructed, x, 1.0e-12_dp, 'reflection modwt/imodwt roundtrip', failures)
   end subroutine test_modwt_roundtrip

   subroutine test_alignment(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions across the test program.

      type(wavelet_transform_type) :: wt
      type(wavelet_transform_type) :: aligned_wt
      type(wavelet_transform_type) :: unaligned_wt
      real(dp) :: x(32, 1)
      real(dp), allocatable :: reconstructed(:, :)
      integer :: ierr
      integer :: i
      integer :: j

      do i = 1, size(x, 1)
         x(i, 1) = real(mod(7 * i, 23), dp) / 8.0_dp
      end do
      call modwt_named(x, wt, filter_name='la8', n_levels=3, ierr=ierr)
      call align(wt, aligned_wt, coe=.false., ierr=ierr)
      call assert_true(ierr == 0 .and. aligned_wt%aligned, 'align status and flag', failures)
      call align(aligned_wt, unaligned_wt, coe=.false., inverse=.true., ierr=ierr)
      call assert_true(ierr == 0 .and. .not. unaligned_wt%aligned, 'inverse align status and flag', failures)
      do j = 1, wt%level
         call assert_close_matrix(unaligned_wt%w(j)%values, wt%w(j)%values, 0.0_dp, 'unaligned W coefficients', failures)
         call assert_close_matrix(unaligned_wt%v(j)%values, wt%v(j)%values, 0.0_dp, 'unaligned V coefficients', failures)
      end do

      call imodwt(aligned_wt, reconstructed, ierr)
      call assert_true(ierr == 0, 'imodwt of aligned object status', failures)
      call assert_close_matrix(reconstructed, x, 1.0e-12_dp, 'aligned imodwt roundtrip', failures)
   end subroutine test_alignment

   subroutine test_mra(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions across the test program.

      type(mra_type) :: result
      real(dp) :: x(64, 1)
      real(dp), allocatable :: reconstructed(:, :)
      integer :: ierr
      integer :: i
      integer :: j

      do i = 1, size(x, 1)
         x(i, 1) = real(mod(11 * i, 29), dp) / 16.0_dp
      end do

      call mra_named(x, result, 3, filter_name='la8', boundary='periodic', method='dwt', ierr=ierr)
      call assert_true(ierr == 0, 'dwt mra status', failures)
      reconstructed = result%s(3)%values
      do j = 1, 3
         reconstructed = reconstructed + result%d(j)%values
      end do
      call assert_close_matrix(reconstructed, x, 2.0e-12_dp, 'dwt mra reconstruction identity', failures)

      call mra_named(x, result, 3, filter_name='la8', boundary='periodic', method='modwt', ierr=ierr)
      call assert_true(ierr == 0, 'modwt mra status', failures)
      reconstructed = result%s(3)%values
      do j = 1, 3
         reconstructed = reconstructed + result%d(j)%values
      end do
      call assert_close_matrix(reconstructed, x, 2.0e-12_dp, 'modwt mra reconstruction identity', failures)
   end subroutine test_mra

   subroutine assert_true(condition, label, failures)
      logical, intent(in) :: condition !! Assertion condition expected to be true.
      character(len=*), intent(in) :: label !! Short diagnostic label identifying the assertion.
      integer, intent(inout) :: failures !! Running failure count incremented when the assertion is false.

      if (.not. condition) then
         failures = failures + 1
         write (*, '(a)') 'FAIL: ' // trim(label)
      end if
   end subroutine assert_true

   subroutine assert_close_vector(actual, expected, tolerance, label, failures)
      real(dp), intent(in) :: actual(:) !! Vector produced by the Fortran implementation.
      real(dp), intent(in) :: expected(:) !! Reference vector expected from the translated operation.
      real(dp), intent(in) :: tolerance !! Maximum permitted absolute elementwise error.
      character(len=*), intent(in) :: label !! Short diagnostic label identifying the assertion.
      integer, intent(inout) :: failures !! Running failure count incremented when the assertion fails.

      if (size(actual) /= size(expected)) then
         failures = failures + 1
         write (*, '(a)') 'FAIL: ' // trim(label) // ' (size mismatch)'
      else if (maxval(abs(actual - expected)) > tolerance) then
         failures = failures + 1
         write (*, '(a,es14.6)') 'FAIL: ' // trim(label) // ' max error=', maxval(abs(actual - expected))
      end if
   end subroutine assert_close_vector

   subroutine assert_close_matrix(actual, expected, tolerance, label, failures)
      real(dp), intent(in) :: actual(:, :) !! Matrix produced by the Fortran implementation.
      real(dp), intent(in) :: expected(:, :) !! Reference matrix expected from the translated operation.
      real(dp), intent(in) :: tolerance !! Maximum permitted absolute elementwise error.
      character(len=*), intent(in) :: label !! Short diagnostic label identifying the assertion.
      integer, intent(inout) :: failures !! Running failure count incremented when the assertion fails.

      if (any(shape(actual) /= shape(expected))) then
         failures = failures + 1
         write (*, '(a)') 'FAIL: ' // trim(label) // ' (shape mismatch)'
      else if (maxval(abs(actual - expected)) > tolerance) then
         failures = failures + 1
         write (*, '(a,es14.6)') 'FAIL: ' // trim(label) // ' max error=', maxval(abs(actual - expected))
      end if
   end subroutine assert_close_matrix

end program test_wavelets
