! SPDX-License-Identifier: GPL-3.0-or-later
program predictive_assessment
   use loo_api
   implicit none
   real(dp) :: ll(400, 8), ll_second(400, 8), response(8), theta(400)
   real(dp) :: scores(8, 2), model_probabilities(2), objective, gap
   real(dp), allocatable :: weights(:, :)
   type(loo_result) :: first, second, waic_result, callback_result
   type(comparison_result) :: comparison
   type(expectation_result) :: expectation
   type(predictive_metric_result) :: metric
   procedure(log_likelihood_callback), pointer :: callback
   integer :: i, j, status
   response = [-1.1_dp, -0.8_dp, -0.2_dp, 0.0_dp, 0.3_dp, 0.5_dp, 0.9_dp, 1.2_dp]
   do i = 1, 400
      ! Deterministic illustrative parameter draws, rather than a fitted MCMC run.
      theta(i) = 0.3_dp * sin(real(i, dp) * 0.37_dp)
      do j = 1, 8
         ll(i, j) = -0.5_dp * (response(j) - theta(i)) ** 2 - 0.9189385332046727_dp
         ll_second(i, j) = -0.5_dp * ((response(j) - theta(i)) / 1.2_dp) ** 2 &
            - log(1.2_dp) - 0.9189385332046727_dp
      end do
   end do
   first = loo(ll)
   second = loo(ll_second)
   waic_result = waic(ll)
   if (first%status /= 0 .or. second%status /= 0 .or. waic_result%status /= 0) error stop 'Invalid example input'
   print '(a,f12.6,a,f12.6)', 'LOO ELPD = ', first%estimates(1, 1), ', SE = ', first%estimates(1, 2)
   print '(a,f12.6)', 'WAIC ELPD = ', waic_result%estimates(1, 1)
   print '(a,8f9.4)', 'Pareto k = ', first%importance%pareto_k
   scores(:, 1) = first%pointwise(:, 1)
   scores(:, 2) = second%pointwise(:, 1)
   comparison = loo_compare(scores)
   print '(a,2i4)', 'Model order = ', comparison%order
   call stacking_weights(scores, model_probabilities, status, objective, kkt_gap = gap)
   if (status /= 0) error stop 'Stacking did not converge'
   print '(a,2f10.6)', 'Stacking weights = ', model_probabilities
   weights = importance_weights(first%importance, logarithm = .false.)
   expectation = E_loo(spread(theta, 2, 8), first%importance, 'mean')
   print '(a,8f9.4)', 'LOO parameter means = ', expectation%value
   metric = loo_predictive_metric(spread(theta, 2, 8), response, weights, 'rmse')
   print '(a,f12.6)', 'Predictive RMSE = ', metric%estimate
   callback => evaluate_log_likelihood
   callback_result = loo_function(callback, 400, 8)
   if (maxval(abs(callback_result%pointwise - first%pointwise)) > 1.0e-12_dp) error stop 'Callback mismatch'
contains
   subroutine evaluate_log_likelihood(observation, values)
      integer, intent(in) :: observation !! One-based observation selected by loo_function.
      real(dp), intent(out) :: values(:) !! Log likelihood for each of the 400 illustrative posterior draws.
      values = ll(:, observation)
   end subroutine evaluate_log_likelihood
end program predictive_assessment
