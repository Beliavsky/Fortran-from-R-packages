program learnbayes_example
   use learnbayes
   implicit none

   type(rng_state) :: rng
   type(blinreg_result) :: fit
   real(dp) :: x(6, 2)
   real(dp) :: y(6)
   real(dp) :: ab(2)
   real(dp) :: predictive(7)
   integer :: info
   integer :: i

   call beta_select(0.05_dp, 0.20_dp, 0.95_dp, 0.60_dp, ab, info)
   if (info /= 0) error stop 'beta prior elicitation failed'
   write (*, '(a,2f10.3)') 'Elicited beta shapes: ', ab

   do i = 0, 6
      predictive(i + 1) = pbetap(ab(1), ab(2), 6, i)
   end do
   write (*, '(a,f10.6)') 'Predictive probabilities sum to: ', sum(predictive)

   x(:, 1) = 1.0_dp
   x(:, 2) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
   y = [1.1_dp, 2.8_dp, 5.2_dp, 6.9_dp, 9.1_dp, 11.2_dp]
   call rng_seed(rng, 20260830_i8)
   call blinreg(rng, y, x, 4000, fit, info)
   if (info /= 0) error stop 'Bayesian linear regression failed'
   write (*, '(a,2f10.4)') 'Posterior mean coefficients: ', &
      sum(fit%beta, dim=1)/real(size(fit%beta, 1), dp)
   write (*, '(a,f10.4)') 'Posterior mean residual sd: ', sum(fit%sigma)/real(size(fit%sigma), dp)
end program learnbayes_example
