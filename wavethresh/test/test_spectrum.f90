! SPDX-License-Identifier: GPL-2.0-or-later
program test_spectrum
   use wavethresh_types, only : dp, wt_vector_t, wd_t, ewspec_result_t
   use wavethresh_spectrum, only : psi_j, psi_j_mat, ipndacw, local_spec_wd, local_spec_wst, ewspec
   implicit none
   real(dp), parameter :: tol = 2.0e-12_dp
   type(wt_vector_t), allocatable :: psi(:)
   real(dp), allocatable :: matrix(:,:)
   real(dp), allocatable :: inner(:,:)
   real(dp) :: expected_haar(4, 4)
   type(wd_t) :: stationary
   type(wd_t) :: periodogram
   type(ewspec_result_t) :: estimate
   real(dp) :: expected_level_zero(4)
   real(dp) :: expected_level_one(4)

   psi = psi_j(-4, 1.0_dp, "DaubExPhase")
   call assert_true(size(psi) == 4, "PsiJ Haar scale count")
   call assert_vector_close(psi(1)%values, [-0.5_dp, 1.0_dp, -0.5_dp], tol, "PsiJ Haar scale 1")
   call assert_vector_close(psi(2)%values, &
      [-0.25_dp, -0.5_dp, 0.25_dp, 1.0_dp, 0.25_dp, -0.5_dp, -0.25_dp], tol, "PsiJ Haar scale 2")
   call assert_true(size(psi(4)%values) == 31, "PsiJ Haar scale 4 length")
   call assert_close(psi(4)%values(16), 1.0_dp, tol, "PsiJ Haar scale 4 center")

   matrix = psi_j_mat(-3, 1.0_dp, "DaubExPhase")
   call assert_true(all(shape(matrix) == [3, 15]), "PsiJmat Haar dimensions")
   call assert_close(matrix(1, 7), -0.5_dp, tol, "PsiJmat Haar first left")
   call assert_close(matrix(1, 8), 1.0_dp, tol, "PsiJmat Haar center")
   call assert_close(matrix(3, 1), -0.125_dp, tol, "PsiJmat Haar third start")
   call assert_close(matrix(3, 15), -0.125_dp, tol, "PsiJmat Haar third end")

   expected_haar(1, :) = [1.5_dp, 0.75_dp, 0.375_dp, 0.1875_dp]
   expected_haar(2, :) = [0.75_dp, 1.75_dp, 1.125_dp, 0.5625_dp]
   expected_haar(3, :) = [0.375_dp, 1.125_dp, 2.875_dp, 2.0625_dp]
   expected_haar(4, :) = [0.1875_dp, 0.5625_dp, 2.0625_dp, 5.4375_dp]
   inner = ipndacw(-4, 1.0_dp, "DaubExPhase")
   call assert_matrix_close(inner, expected_haar, tol, "ipndacw Haar")

   psi = psi_j(-1, 10.0_dp, "DaubLeAsymm")
   call assert_true(size(psi) == 1 .and. size(psi(1)%values) == 39, "PsiJ LA10 scale 1 dimensions")
   call assert_close(psi(1)%values(1), 3.537571e-7_dp, 5.0e-13_dp, "PsiJ LA10 first coefficient")
   call assert_close(psi(1)%values(19), -6.209080e-1_dp, 5.0e-8_dp, "PsiJ LA10 near center")
   call assert_close(psi(1)%values(20), 1.0_dp, 1.0e-9_dp, "PsiJ LA10 center")

   inner = ipndacw(-2, 10.0_dp, "DaubLeAsymm")
   call assert_close(inner(1, 1), 1.839101_dp, 1.0e-6_dp, "ipndacw LA10 11")
   call assert_close(inner(1, 2), 0.3215934_dp, 1.0e-6_dp, "ipndacw LA10 12")
   call assert_close(inner(2, 2), 3.035353_dp, 2.0e-6_dp, "ipndacw LA10 22")

   stationary%n_original = 4
   stationary%nlevels = 2
   stationary%transform_type = "station"
   stationary%ok = .true.
   stationary%message = "ok"
   allocate(stationary%detail(0:1))
   stationary%detail(0)%values = [-2.0_dp, -1.0_dp, 0.5_dp, 3.0_dp]
   stationary%detail(1)%values = [4.0_dp, -3.0_dp, 2.0_dp, -1.0_dp]
   periodogram = local_spec_wd(stationary)
   call assert_true(periodogram%ok, "LocalSpec.wd raw periodogram status")
   call assert_vector_close(periodogram%detail(0)%values, [4.0_dp, 1.0_dp, 0.25_dp, 9.0_dp], tol, &
      "LocalSpec.wd level zero")
   call assert_vector_close(periodogram%detail(1)%values, [16.0_dp, 9.0_dp, 4.0_dp, 1.0_dp], tol, &
      "LocalSpec.wd level one")
   periodogram = local_spec_wst(stationary)
   call assert_vector_close(periodogram%detail(0)%values, [4.0_dp, 1.0_dp, 0.25_dp, 9.0_dp], tol, &
      "LocalSpec.wst dispatch")

   estimate = ewspec(stationary, 1.0_dp, "DaubExPhase")
   call assert_true(estimate%ok, "ewspec unsmoothed stationary status")
   expected_level_one = (28.0_dp * [16.0_dp, 9.0_dp, 4.0_dp, 1.0_dp] - &
      12.0_dp * [4.0_dp, 1.0_dp, 0.25_dp, 9.0_dp]) / 33.0_dp
   expected_level_zero = (-12.0_dp * [16.0_dp, 9.0_dp, 4.0_dp, 1.0_dp] + &
      24.0_dp * [4.0_dp, 1.0_dp, 0.25_dp, 9.0_dp]) / 33.0_dp
   call assert_vector_close(estimate%spectrum%detail(1)%values, expected_level_one, tol, &
      "ewspec corrected fine level")
   call assert_vector_close(estimate%spectrum%detail(0)%values, expected_level_zero, tol, &
      "ewspec corrected coarse level")
   call assert_matrix_close(estimate%rm, reshape([1.5_dp, 0.75_dp, 0.75_dp, 1.75_dp], [2, 2]), tol, &
      "ewspec Haar correction matrix")

   print *, "test_spectrum: PASS"

contains

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition !! Condition that must evaluate true for the test to pass.
      character(len=*), intent(in) :: message !! Diagnostic label printed when the assertion fails.
      if (.not. condition) then
         print *, "FAIL: ", trim(message)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual !! Scalar value produced by the translated implementation.
      real(dp), intent(in) :: expected !! Deterministic scalar reference value.
      real(dp), intent(in) :: tolerance !! Maximum permitted absolute difference from the reference value.
      character(len=*), intent(in) :: message !! Diagnostic label printed when the assertion fails.
      if (abs(actual - expected) > tolerance) then
         print *, "FAIL: ", trim(message), actual, expected
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_vector_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual(:) !! Vector produced by the translated implementation.
      real(dp), intent(in) :: expected(:) !! Deterministic reference vector.
      real(dp), intent(in) :: tolerance !! Maximum permitted elementwise absolute difference.
      character(len=*), intent(in) :: message !! Diagnostic label printed when the assertion fails.
      call assert_true(size(actual) == size(expected), trim(message)//" size")
      if (size(actual) > 0) call assert_close(maxval(abs(actual - expected)), 0.0_dp, tolerance, message)
   end subroutine assert_vector_close

   subroutine assert_matrix_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual(:,:) !! Matrix produced by the translated implementation.
      real(dp), intent(in) :: expected(:,:) !! Deterministic reference matrix.
      real(dp), intent(in) :: tolerance !! Maximum permitted elementwise absolute difference.
      character(len=*), intent(in) :: message !! Diagnostic label printed when the assertion fails.
      call assert_true(all(shape(actual) == shape(expected)), trim(message)//" shape")
      if (size(actual) > 0) call assert_close(maxval(abs(actual - expected)), 0.0_dp, tolerance, message)
   end subroutine assert_matrix_close

end program test_spectrum
