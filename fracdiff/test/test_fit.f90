! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran port of fracdiff; see NOTICE.md for attribution.

program test_fit
   use, intrinsic :: iso_fortran_env, only : int64
   use fracdiff, only : dp, fd_ok, fd_iteration_limit, fracdiff_simulation, fracdiff_model, &
                        fracdiff_summary, fracdiff_sim, fracdiff_fit, fracdiff_var, &
                        fracdiff_coefficients, fracdiff_confint, summarize_fracdiff, &
                        fracdiff_aic, fracdiff_bic
   use test_support
   implicit none

   type(fracdiff_simulation) :: simulated
   type(fracdiff_model) :: model, model0
   type(fracdiff_summary) :: summary
   real(dp),allocatable :: coefficients(:),intervals(:,:)
   real(dp) :: ar(1),ma(1),old_h
   integer :: status

   ar=0.2_dp
   ma=-0.4_dp
   simulated=fracdiff_sim(1200,0.3_dp,ar=ar,ma=ma,seed=107_int64,n_start=100)
   model=fracdiff_fit(simulated%series,nar=1,nma=1)
   call assert_true(model%status==fd_ok .or. model%status==fd_iteration_limit,"ARFIMA fit status")
   call assert_close(model%d,0.3_dp,0.12_dp,"recover d")
   call assert_close(model%ar(1),0.2_dp,0.16_dp,"recover AR")
   call assert_close(model%ma(1),-0.4_dp,0.14_dp,"recover MA")
   call assert_close(model%sigma,1.0_dp,0.08_dp,"recover innovation sigma")
   call assert_matrix_symmetric(model%hessian,1.0e-10_dp,"Hessian symmetry")
   call assert_matrix_symmetric(model%covariance,1.0e-10_dp,"covariance symmetry")
   call assert_true(all(model%std_error>0.0_dp),"positive standard errors")

   coefficients=fracdiff_coefficients(model)
   call assert_true(size(coefficients)==3,"coefficient count")
   intervals=fracdiff_confint(model)
   call assert_true(size(intervals,1)==3 .and. size(intervals,2)==2,"confidence interval dimensions")
   call assert_true(all(intervals(:,1)<intervals(:,2)),"ordered confidence intervals")
   summary=summarize_fracdiff(model)
   call assert_close(summary%aic,fracdiff_aic(model),1.0e-14_dp,"summary AIC")
   call assert_close(summary%bic,fracdiff_bic(model),1.0e-14_dp,"summary BIC")

   old_h=model%h
   call fracdiff_var(simulated%series,model,old_h*0.75_dp,status)
   call assert_true(status==fd_ok,"variance recomputation")
   call assert_close(model%h,old_h*0.75_dp,1.0e-15_dp,"updated finite-difference step")

   simulated=fracdiff_sim(1000,0.25_dp,seed=8_int64)
   model0=fracdiff_fit(simulated%series)
   call assert_true(model0%status==fd_ok,"fractional-noise fit status")
   call assert_close(model0%d,0.25_dp,0.10_dp,"fractional-noise d recovery")

   write(*,'(a)') "test_fit: PASS"
end program test_fit
