! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran port of fracdiff; see NOTICE.md for attribution.

program fracdiff_demo
   use, intrinsic :: iso_fortran_env, only : int64
   use fracdiff, only : dp, fracdiff_simulation, fracdiff_model, fracdiff_summary, &
                        fracdiff_sim, fracdiff_fit, summarize_fracdiff
   implicit none

   type(fracdiff_simulation) :: simulated
   type(fracdiff_model) :: fitted
   type(fracdiff_summary) :: summary
   real(dp) :: ar(1), ma(1)
   integer :: i

   ar = 0.2_dp
   ma = -0.4_dp
   simulated = fracdiff_sim(1200, 0.3_dp, ar=ar, ma=ma, seed=107_int64, n_start=100)
   if (size(simulated%series) == 0) error stop "simulation failed"

   fitted = fracdiff_fit(simulated%series, nar=1, nma=1)
   summary = summarize_fracdiff(fitted)

   write(*,'(a)') "ARFIMA(1,d,1) demonstration"
   write(*,'(a,f12.6)') "log likelihood: ", fitted%log_likelihood
   write(*,'(a,f12.6)') "sigma:          ", fitted%sigma
   write(*,'(a)') "parameter       estimate      std. error"
   do i = 1, size(summary%coefficients,1)
      write(*,'(i5,2f15.7)') i, summary%coefficients(i,1), summary%coefficients(i,2)
   end do
   write(*,'(a,f12.4)') "AIC: ", summary%aic
   write(*,'(a,a)') "status: ", fitted%message
end program fracdiff_demo
