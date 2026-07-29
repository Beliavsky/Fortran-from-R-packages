! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran port of fracdiff; see NOTICE.md for attribution.

program test_simulation
   use, intrinsic :: iso_fortran_env, only : int64
   use fracdiff, only : dp, fd_ok, fracdiff_simulation, fracdiff_sim, &
                        fractional_arma_filter_simulation
   use test_support
   implicit none

   real(dp) :: innovations(13), ar(1), ma(1), series(12), expected(12), shifted(12)
   type(fracdiff_simulation) :: simulated_a, simulated_b
   integer :: i, status

   do i=1,size(innovations)
      innovations(i)=sin(real(i,dp))
   end do
   ar=0.5_dp
   ma=-0.25_dp
   expected=[1.3638514015689536_dp,1.4621770995202188_dp,0.37115787477770756_dp, &
      -0.8406744038598260_dp,-1.0434047084106381_dp,-0.05823951658320525_dp, &
      1.1926778790566490_dp,1.5420015683652310_dp,0.6538108897322745_dp, &
      -0.6671469321103923_dp,-1.2162988171678542_dp,-0.4977034349143774_dp]
   call fractional_arma_filter_simulation(innovations,ar,ma,0.2_dp,0.0_dp,series,status)
   call assert_true(status==fd_ok,"fixed-innovation simulation status")
   call assert_vector_close(series,expected,2.0e-12_dp,"fixed-innovation simulation reference")
   call fractional_arma_filter_simulation(innovations,ar,ma,0.2_dp,2.0_dp,shifted,status)
   call assert_vector_close(shifted,expected+2.0_dp,2.0e-12_dp,"simulation mean added after filtering")

   simulated_a=fracdiff_sim(200,0.25_dp,ar=ar,ma=ma,seed=123_int64,n_start=20)
   simulated_b=fracdiff_sim(200,0.25_dp,ar=ar,ma=ma,seed=123_int64,n_start=20)
   call assert_true(simulated_a%status==fd_ok,"seeded simulation status")
   call assert_vector_close(simulated_a%series,simulated_b%series,0.0_dp,"deterministic seed")
   call assert_true(maxval(abs(simulated_a%series))<20.0_dp,"simulated values finite and plausible")

   write(*,'(a)') "test_simulation: PASS"
end program test_simulation
