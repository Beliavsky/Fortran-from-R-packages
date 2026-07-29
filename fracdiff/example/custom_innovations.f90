! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran port of fracdiff; see NOTICE.md for attribution.

program custom_innovations
   use fracdiff, only : dp, fracdiff_simulation, fracdiff_sim
   implicit none

   integer, parameter :: n = 20, n_start = 4
   real(dp) :: innovations(n+1), start_innovations(n_start)
   real(dp) :: ar(1), ma(1)
   type(fracdiff_simulation) :: simulated
   integer :: i

   ar = 0.5_dp
   ma = -0.25_dp
   do i = 1, size(innovations)
      innovations(i) = sin(real(i,dp))
   end do
   do i = 1, n_start
      start_innovations(i) = cos(real(i,dp))
   end do

   simulated = fracdiff_sim(n, 0.2_dp, ar=ar, ma=ma, n_start=n_start, &
      innovations=innovations, start_innovations=start_innovations)
   write(*,'(a)') "First ten simulated observations:"
   write(*,'(5f14.7)') simulated%series(1:10)
end program custom_innovations
