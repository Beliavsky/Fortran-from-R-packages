! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran port of fracdiff; see NOTICE.md for attribution.

program test_estimators
   use fracdiff, only : dp, fd_ok, fractional_d_estimate, fd_gph, fd_sperio
   use test_support
   implicit none

   integer, parameter :: n=256
   real(dp) :: x(n)
   type(fractional_d_estimate) :: gph, sperio
   integer :: i

   do i=1,n
      x(i)=sin(0.17_dp*real(i,dp))+0.35_dp*cos(0.031_dp*real(i,dp))+0.002_dp*real(i,dp)
   end do
   gph=fd_gph(x)
   sperio=fd_sperio(x)
   call assert_true(gph%status==fd_ok,"GPH status")
   call assert_true(sperio%status==fd_ok,"Sperio status")
   call assert_close(gph%d,1.0078245382412603_dp,2.0e-12_dp,"GPH d reference")
   call assert_close(gph%sd_asymptotic,0.2102807539782185_dp,2.0e-12_dp,"GPH asymptotic SE")
   call assert_close(gph%sd_regression,0.34897788608158936_dp,2.0e-12_dp,"GPH regression SE")
   call assert_close(sperio%d,0.9617721093864684_dp,2.0e-12_dp,"Sperio d reference")
   call assert_close(sperio%sd_asymptotic,0.09123744883801196_dp,2.0e-12_dp,"Sperio asymptotic SE")
   call assert_close(sperio%sd_regression,0.30917577437627675_dp,2.0e-12_dp,"Sperio regression SE")

   write(*,'(a)') "test_estimators: PASS"
end program test_estimators
