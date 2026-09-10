! SPDX-License-Identifier: GPL-2.0-or-later
! Test/example code for the modern Fortran translation of the upstream Rssa package.
program test_iterative_oblique
   use rssa, only : dp, eossa, iossa, reconstruct_ssa, rssa_success, ssa_1d, ssa_result
   implicit none
   real(dp) :: expected_first(12), expected_second(12), source(12)
   real(dp), allocatable :: component(:), reconstruction(:)
   type(ssa_result) :: eos, fit, ios
   integer, allocatable :: clusters(:)
   integer :: i, iterations
   logical :: converged

   do i = 1, size(source)
      expected_first(i) = 2.0_dp * 0.82_dp ** real(i - 1, dp)
      expected_second(i) = 0.7_dp * (-0.45_dp) ** real(i - 1, dp)
   end do
   source = expected_first + expected_second
   fit = ssa_1d(source, window=6, neig=2)
   call check(fit%info == rssa_success, "ordinary rank-2 SSA setup")

   ios = iossa(fit, indices=[1, 2], group_labels=[1, 2], tol=1.0e-10_dp, maxiter=50, &
               iterations=iterations, converged=converged)
   call check(ios%info == rssa_success, "I-OSSA status")
   call check(converged, "I-OSSA convergence")
   call check(iterations <= 50, "I-OSSA iteration bound")
   component = reconstruct_ssa(ios, [1])
   call check(maxval(abs(component - expected_first)) < 2.0e-7_dp, "I-OSSA first exponential")
   component = reconstruct_ssa(ios, [2])
   call check(maxval(abs(component - expected_second)) < 2.0e-7_dp, "I-OSSA second exponential")
   reconstruction = reconstruct_ssa(ios, [1, 2])
   call check(maxval(abs(reconstruction - source)) < 2.0e-10_dp, "I-OSSA trajectory preservation")

   eos = eossa(fit, indices=[1, 2], cluster_labels=clusters)
   if (eos%info /= rssa_success) print '(a,i0)', "EOSSA status code: ", eos%info
   call check(eos%info == rssa_success, "EOSSA status")
   call check(size(clusters) == 2, "EOSSA cluster count")
   call check(all(clusters == [1, 2]), "EOSSA frequency-ordered clusters")
   component = reconstruct_ssa(eos, [1])
   call check(maxval(abs(component - expected_first)) < 2.0e-7_dp, "EOSSA first exponential")
   component = reconstruct_ssa(eos, [2])
   call check(maxval(abs(component - expected_second)) < 2.0e-7_dp, "EOSSA second exponential")
   reconstruction = reconstruct_ssa(eos, [1, 2])
   call check(maxval(abs(reconstruction - source)) < 2.0e-10_dp, "EOSSA trajectory preservation")

   print '(a)', 'test_iterative_oblique: PASS'
contains
   subroutine check(condition, message)
      logical, intent(in) :: condition !! Assertion predicate that must be true.
      character(len=*), intent(in) :: message !! Short description printed when the assertion fails.
      if (.not. condition) then
         print '(a)', 'test_iterative_oblique: FAIL: ' // trim(message)
         error stop 1
      end if
   end subroutine check
end program test_iterative_oblique
