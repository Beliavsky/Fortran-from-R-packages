! SPDX-License-Identifier: GPL-2.0-or-later
program test_pairs_location
   use icsnp, only : dp, icsnp_ok, pair_diff, pair_sum, pair_prod, spatial_median, &
      hl_loc, vdw_loc
   implicit none
   real(dp) :: x(5, 2), v(5), loc
   real(dp), allocatable :: a(:,:), center(:)
   integer :: status, iterations

   x = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, &
                2.0_dp, 1.0_dp, 4.0_dp, 3.0_dp, 5.0_dp], shape(x))
   call pair_diff(x, a, status)
   call check(status == icsnp_ok .and. size(a, 1) == 10, 'pair_diff size')
   call check(maxval(abs(a(1, :) - (x(1, :) - x(2, :)))) < 1.0e-12_dp, 'pair_diff values')
   call pair_sum(x, a, status)
   call check(maxval(abs(a(1, :) - (x(1, :) + x(2, :)))) < 1.0e-12_dp, 'pair_sum values')
   call pair_prod(x, a, status)
   call check(maxval(abs(a(1, :) - x(1, :) * x(2, :))) < 1.0e-12_dp, 'pair_prod values')

   call spatial_median(x, center, status, iterations)
   call check(status == icsnp_ok, 'spatial_median status')
   call check(size(center) == 2 .and. all(abs(center - [3.0_dp, 3.0_dp]) < 0.8_dp), &
      'spatial_median value')

   v = [1.0_dp, 2.0_dp, 3.0_dp, 9.0_dp, 10.0_dp]
   loc = hl_loc(v, status)
   call check(status == icsnp_ok .and. loc > 2.5_dp .and. loc < 7.0_dp, 'hl_loc')
   loc = vdw_loc(v, status)
   call check(status == icsnp_ok .and. loc > minval(v) .and. loc < maxval(v), 'vdw_loc')
   print '(a)', 'test_pairs_location: PASS'
contains
   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a)') 'FAIL: ' // trim(label)
         error stop 1
      end if
   end subroutine check
end program test_pairs_location
