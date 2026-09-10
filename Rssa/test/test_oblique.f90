! SPDX-License-Identifier: GPL-2.0-or-later
! Test/example code for the modern Fortran translation of the upstream Rssa package.
program test_oblique
   use rssa_kinds, only : dp
   use rssa_oblique, only : decompose_wossa, owcor_ssa, wcor_ossa
   use rssa_types, only : rssa_success, ssa_result
   implicit none
   type(ssa_result) :: object, weighted
   real(dp), allocatable :: matrix(:, :), ordinary(:, :)
   real(dp), parameter :: root2 = sqrt(2.0_dp)
   real(dp), parameter :: root3 = sqrt(3.0_dp)
   real(dp), parameter :: root6 = sqrt(6.0_dp)
   real(dp) :: expected_u(2), expected_v(4), series(5)
   integer :: info

   series = [1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, 16.0_dp]
   weighted = decompose_wossa(series, 2, neig=1, &
                              column_weights=[4.0_dp, 0.25_dp], &
                              row_weights=[1.0_dp, 4.0_dp, 0.25_dp, 9.0_dp])
   if (weighted%info /= rssa_success) then
      print '(a,i0)', "decompose_wossa status: ", weighted%info
      error stop "decompose_wossa status"
   end if
   if (size(weighted%sigma) /= 1) error stop "decompose_wossa rank"
   if (abs(weighted%sigma(1) - sqrt(425.0_dp)) > 2.0e-12_dp) error stop "decompose_wossa sigma"
   expected_u = [1.0_dp, 2.0_dp] / sqrt(5.0_dp)
   expected_v = [1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp] / sqrt(85.0_dp)
   call assert_close(abs(weighted%u(:, 1)), expected_u, 2.0e-12_dp, "decompose_wossa U")
   call assert_close(abs(weighted%v(:, 1)), expected_v, 2.0e-12_dp, "decompose_wossa V")

   allocate(object%series(4), object%sigma(2), object%u(2, 2), object%v(3, 2))
   object%series = 0.0_dp
   object%window = 2
   object%sigma = [3.0_dp, 1.5_dp]
   object%u(:, 1) = [root3 / 2.0_dp, 0.5_dp]
   object%u(:, 2) = [-0.5_dp, root3 / 2.0_dp]
   object%v(:, 1) = [1.0_dp / root2, 1.0_dp / root2, 0.0_dp]
   object%v(:, 2) = [1.0_dp / root6, -1.0_dp / root6, 2.0_dp / root6]

   matrix = owcor_ssa(object, info=info)
   if (info /= rssa_success) error stop "owcor_ssa status"
   if (maxval(abs(matrix - reshape([1.0_dp, 0.0139234562739617_dp, &
                                    0.0139234562739617_dp, 1.0_dp], [2, 2]))) > 3.0e-13_dp) then
      error stop "owcor_ssa reference"
   end if

   ordinary = wcor_ossa(object, info=info)
   if (info /= rssa_success) error stop "wcor_ossa status"
   if (maxval(abs(ordinary - reshape([1.0_dp, 0.0102341742038767_dp, &
                                      0.0102341742038767_dp, 1.0_dp], [2, 2]))) > 3.0e-13_dp) then
      error stop "wcor_ossa reference"
   end if

   matrix = owcor_ssa(object, group_labels=[1, 1], info=info)
   if (info /= rssa_success) error stop "owcor_ssa grouped status"
   if (any(shape(matrix) /= [1, 1])) error stop "owcor_ssa grouped shape"
   if (abs(matrix(1, 1) - 1.0_dp) > 2.0e-14_dp) error stop "owcor_ssa grouped value"

   print *, "test_oblique: PASS"
contains
   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual(:) !! Values produced by the translated implementation.
      real(dp), intent(in) :: expected(:) !! Independently calculated reference values.
      real(dp), intent(in) :: tolerance !! Maximum permitted absolute error.
      character(len=*), intent(in) :: label !! Failure label identifying the checked quantity.

      if (size(actual) /= size(expected)) error stop label // ": size mismatch"
      if (size(actual) > 0) then
         if (maxval(abs(actual - expected)) > tolerance) error stop label // ": value mismatch"
      end if
   end subroutine assert_close
end program test_oblique
