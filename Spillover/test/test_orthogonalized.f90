! SPDX-License-Identifier: GPL-2.0-only
! Modern Fortran translation of computational code from the Spillover R package.
program test_orthogonalized
   use spillover
   implicit none

   real(dp) :: ar(2, 2, 1), sigma(2, 2)
   real(dp), allocatable :: fevd(:, :, :), table(:, :)
   type(var_model) :: model
   type(orthogonal_average_result) :: average
   type(spillover_result) :: result, source_result
   integer :: info
   character(len=200) :: message

   ar = 0.0_dp
   sigma = reshape([1.0_dp, 0.5_dp, 0.5_dp, 1.0_dp], [2, 2])
   call initialize_var_model(ar, sigma, model, info=info, message=message)
   call assert_true(info == spillover_success, trim(message))

   call orthogonalized_fevd(model, 1, fevd, info=info, message=message)
   call assert_true(info == spillover_success, trim(message))
   call assert_close(fevd(1, 1, 1), 1.0_dp, 1.0e-13_dp, 'single own share first variable')
   call assert_close(fevd(1, 2, 1), 0.0_dp, 1.0e-13_dp, 'single cross share first variable')
   call assert_close(fevd(2, 1, 1), 0.25_dp, 1.0e-13_dp, 'single cross share second variable')
   call assert_close(fevd(2, 2, 1), 0.75_dp, 1.0e-13_dp, 'single own share second variable')

   call orthogonal_average_exact(model, 1, average, info=info, message=message)
   call assert_true(info == spillover_success, trim(message))
   call assert_true(average%n_permutations == 2, 'two exact permutations')
   call assert_close(average%average(1, 1), 87.5_dp, 1.0e-12_dp, 'average diagonal')
   call assert_close(average%average(1, 2), 12.5_dp, 1.0e-12_dp, 'average off diagonal')

   call orthogonalized_spillover(model, 1, ortho_total, .true., result, info=info, message=message)
   call assert_true(info == spillover_success, trim(message))
   call assert_close(result%total, 12.5_dp, 1.0e-12_dp, 'exact-order total index')

   call orthogonalized_spillover(model, 1, ortho_partial, .true., source_result, &
      source_compatible=.true., info=info, message=message)
   call assert_true(info == spillover_success, trim(message))
   call assert_close(maxval(abs(source_result%shares - result%shares)), 0.0_dp, 1.0e-12_dp, &
      'source-compatible partial/total label swap')

   call compatibility_table(result, table, info, message)
   call assert_true(info == spillover_success, trim(message))
   call assert_true(size(table, 1) == 4 .and. size(table, 2) == 3, 'compatibility table shape')
   call assert_close(table(3, 3), result%total, 1.0e-12_dp, 'compatibility total cell')
   call assert_close(table(4, 3), 100.0_dp, 1.0e-12_dp, 'compatibility grand total')

   print '(a)', 'test_orthogonalized: PASS'

contains

   subroutine assert_true(condition, text)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: text
      if (.not. condition) then
         write(*, '(a)') 'FAIL: ' // text
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(actual, expected, tolerance, text)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: text
      if (abs(actual - expected) > tolerance) then
         write(*, '(a,2(1x,es14.6))') 'FAIL: ' // text, actual, expected
         error stop 1
      end if
   end subroutine assert_close

end program test_orthogonalized
