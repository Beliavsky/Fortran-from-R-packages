! SPDX-License-Identifier: GPL-2.0-or-later
program test_transforms
   use moments, only : dp, all_moments, raw2central, central2raw, all_cumulants
   implicit none

   real(dp) :: x(4), matrix_raw(5, 2), eta(2)
   real(dp), allocatable :: raw(:), central(:), recovered(:), cumulants(:), legacy(:)
   real(dp), allocatable :: matrix_central(:, :), matrix_recovered(:, :), matrix_cumulants(:, :)

   x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
   raw = all_moments(x, 4)
   central = raw2central(raw)
   call check_close(central(1), 1.0_dp, 1.0e-14_dp, 'central order zero')
   call check_close(central(2), 0.0_dp, 1.0e-14_dp, 'central order one')
   call check_close(central(3), 1.25_dp, 1.0e-14_dp, 'central order two')
   call check_close(central(4), 0.0_dp, 1.0e-14_dp, 'central order three')
   call check_close(central(5), 2.5625_dp, 1.0e-14_dp, 'central order four')

   recovered = central2raw(central, raw(2))
   call check(maxval(abs(recovered - raw)) < 1.0e-13_dp, 'central/raw round trip')

   cumulants = all_cumulants(raw)
   call check_close(cumulants(1), 0.0_dp, 1.0e-14_dp, 'cumulant zero')
   call check_close(cumulants(2), 2.5_dp, 1.0e-14_dp, 'cumulant one')
   call check_close(cumulants(3), 1.25_dp, 1.0e-14_dp, 'cumulant two')
   call check_close(cumulants(4), 0.0_dp, 1.0e-13_dp, 'cumulant three')
   call check_close(cumulants(5), -2.125_dp, 1.0e-12_dp, 'cumulant four')

   legacy = all_cumulants(raw, legacy=.true.)
   call check_close(legacy(2), 0.0_dp, 1.0e-14_dp, 'legacy first cumulant')
   call check_close(legacy(4), 18.75_dp, 1.0e-12_dp, 'legacy third cumulant')

   matrix_raw(:, 1) = raw
   matrix_raw(:, 2) = [1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, 16.0_dp] * raw
   eta = [2.5_dp, 5.0_dp]
   matrix_central = raw2central(matrix_raw)
   matrix_recovered = central2raw(matrix_central, eta)
   call check(maxval(abs(matrix_recovered - matrix_raw)) < 1.0e-11_dp, 'matrix round trip')
   matrix_cumulants = all_cumulants(matrix_raw)
   call check_close(matrix_cumulants(3, 2), 5.0_dp, 1.0e-12_dp, 'matrix cumulants')

   print '(a)', 'test_transforms: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a,1x,a)') 'FAIL:', trim(label)
         error stop 1
      end if
   end subroutine check

   subroutine check_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label
      call check(abs(actual - expected) <= tolerance, label)
   end subroutine check_close

end program test_transforms
