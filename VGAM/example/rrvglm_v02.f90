program rrvglm_v02
   use vgam
   implicit none
   integer, parameter :: n = 40, m = 3
   real(dp) :: x(n, 3), y(n, m), t
   integer :: i
   integer :: families(m)
   type(rrvglm_result_t) :: fit

   do i = 1, n
      t = -1.0_dp + 2.0_dp*real(i - 1, dp)/real(n - 1, dp)
      x(i, :) = [1.0_dp, t, sin(2.0_dp*t)]
   end do
   y(:, 1) = 0.3_dp + 0.8_dp*x(:, 2) - 0.4_dp*x(:, 3)
   y(:, 2) = -0.2_dp - 0.4_dp*x(:, 2) + 0.2_dp*x(:, 3)
   y(:, 3) = 0.5_dp + 1.2_dp*x(:, 2) - 0.6_dp*x(:, 3)
   families = family_gaussian

   call fit_rrvglm(y, x, 1, families, fit)
   print '(a,l1)', 'RR-VGLM converged: ', fit%converged
   print '(a,i0)', 'Reduced coefficient rank: ', fit%effective_rank
   print '(a,f12.6)', 'Residual deviance: ', fit%deviance
end program rrvglm_v02
