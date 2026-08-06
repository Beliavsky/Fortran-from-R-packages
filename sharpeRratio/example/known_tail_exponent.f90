! SPDX-License-Identifier: GPL-3.0-only
program known_tail_exponent
   use, intrinsic :: iso_fortran_env, only : int64
   use sharpe_rratio, only : dp, snr_result, estimate_snr
   implicit none
   real(dp) :: returns(100)
   type(snr_result) :: estimate
   integer :: i

   do i = 1, size(returns)
      returns(i) = 0.002_dp+0.01_dp*sin(0.37_dp*real(i,dp))
   end do
   estimate = estimate_snr(returns,num_perm=30,nu=5.0_dp,seed=17_int64)
   if (.not. estimate%ok) error stop trim(estimate%message)
   print '(a,f12.6)', 'moment-free SNR: ',estimate%snr
end program known_tail_exponent
