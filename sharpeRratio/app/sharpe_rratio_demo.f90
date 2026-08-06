! SPDX-License-Identifier: GPL-3.0-only
program sharpe_rratio_demo
   use, intrinsic :: iso_fortran_env, only : int64
   use ghyp, only : ghyp_model_type, student_t_uv, rghyp
   use ghyp_kinds, only : i8
   use sharpe_rratio, only : dp, snr_result, estimate_snr
   implicit none
   type(ghyp_model_type) :: model
   type(snr_result) :: estimate
   real(dp), allocatable :: returns(:,:)
   logical :: ok

   model = student_t_uv(3.0_dp,mu=0.05_dp)
   call rghyp(250,model,returns,ok,20260803_i8)
   if (.not. ok) error stop 'simulation failed'

   estimate = estimate_snr(returns(:,1),num_perm=50,seed=42_int64, &
      max_fit_iterations=500)
   if (.not. estimate%ok) error stop trim(estimate%message)

   print '(a,f12.6)', 'SNR             : ',estimate%snr
   print '(a,2f12.6)', 'confidence range: ',estimate%ci_lower,estimate%ci_upper
   print '(a,f12.4)', 'tail exponent   : ',estimate%nu
   print '(a,f12.4)', 'mean R0         : ',estimate%r0bar
   print '(a,i0)', 'observations    : ',estimate%n
end program sharpe_rratio_demo
