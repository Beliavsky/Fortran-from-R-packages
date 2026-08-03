! SPDX-License-Identifier: GPL-2.0-or-later
program location_example
   use icsnp, only : dp, icsnp_ok, spatial_median, hl_loc, vdw_loc
   implicit none
   real(dp) :: x(9, 2), values(9)
   real(dp), allocatable :: center(:)
   real(dp) :: hl, vdw
   integer :: i, status, iterations

   do i = 1, 9
      x(i, 1) = real(i - 5, dp)
      x(i, 2) = 0.4_dp * real(i - 5, dp) + sin(real(i, dp))
   end do
   x(9, :) = [15.0_dp, -10.0_dp]
   values = x(:, 1)

   call spatial_median(x, center, status, iterations)
   if (status /= icsnp_ok) error stop 'spatial_median failed'
   hl = hl_loc(values, status)
   if (status /= icsnp_ok) error stop 'hl_loc failed'
   vdw = vdw_loc(values, status)
   if (status /= icsnp_ok) error stop 'vdw_loc failed'

   write(*, '(a,2f12.6)') 'Spatial median: ', center
   write(*, '(a,f12.6)') 'Hodges-Lehmann location: ', hl
   write(*, '(a,f12.6)') 'Van der Waerden location: ', vdw
end program location_example
