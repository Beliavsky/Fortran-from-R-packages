! SPDX-License-Identifier: GPL-2.0-or-later
! Test/example code for the modern Fortran translation of the upstream Rssa package.
program test_decomposition
   use rssa, only : cssa_result, decompose_pssa, dp, mssa_result, reconstruct_complex, reconstruct_mssa
   use rssa, only : reconstruct_ssa, ssa_1d, ssa_complex, ssa_mssa, ssa_result
   implicit none
   real(dp) :: x(8), channels(8, 2), projector(4, 1)
   real(dp), allocatable :: reconstructed(:), multi(:, :)
   complex(dp) :: z(8)
   complex(dp), allocatable :: z_reconstructed(:)
   type(ssa_result) :: fit, pfit
   type(mssa_result) :: mfit
   type(cssa_result) :: zfit
   integer :: i, info
   x = [(2.0_dp**real(i - 1, dp), i=1, 8)]
   fit = ssa_1d(x, 4, 1)
   reconstructed = reconstruct_ssa(fit, [1])
   call assert_true(fit%info == 0, 'SSA status')
   call assert_close(maxval(abs(reconstructed - x)), 0.0_dp, 1.0e-8_dp, 'rank-one reconstruction')
   channels(:, 1) = x
   channels(:, 2) = 3.0_dp * x
   mfit = ssa_mssa(channels, 4, 1)
   multi = reconstruct_mssa(mfit, [1])
   call assert_true(mfit%info == 0 .and. all(shape(multi) == shape(channels)), 'MSSA reconstruction')
   projector = 0.5_dp
   pfit = decompose_pssa(x, 4, column_projector=projector, neig=1, info=info)
   call assert_true(info == 0 .and. allocated(pfit%sigma), 'projection SSA')
   do i = 1, size(z)
      z(i) = cmplx(0.8_dp, 0.3_dp, dp)**(i - 1)
   end do
   zfit = ssa_complex(z, 4, 1)
   z_reconstructed = reconstruct_complex(zfit, [1])
   call assert_true(zfit%info == 0, 'complex SSA status')
   call assert_close(maxval(abs(z_reconstructed - z)), 0.0_dp, 1.0e-10_dp, 'complex rank-one reconstruction')
   print *, 'test_decomposition: PASS'
contains
   subroutine assert_true(condition, message)
      logical, intent(in) :: condition !! Condition required to be true.
      character(len=*), intent(in) :: message !! Failure diagnostic.
      if (.not. condition) error stop message
   end subroutine assert_true
   subroutine assert_close(value, expected, tolerance, message)
      real(dp), intent(in) :: value !! Computed scalar value.
      real(dp), intent(in) :: expected !! Expected scalar reference.
      real(dp), intent(in) :: tolerance !! Maximum absolute error.
      character(len=*), intent(in) :: message !! Failure diagnostic.
      if (abs(value - expected) > tolerance) error stop message
   end subroutine assert_close
end program test_decomposition
