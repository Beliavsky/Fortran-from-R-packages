! SPDX-License-Identifier: GPL-3.0-or-later
program test_core
   use loo_api
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan, &
      ieee_positive_inf, ieee_negative_inf
   implicit none
   real(dp) :: ll(200, 4), eff(4), lp(200), lg(200), x(200, 4), a(50, 4, 4)
   real(dp) :: shifted(200, 4), c(4, 3), q(3), t(100), values(1), nan, inf
   real(dp), allocatable :: w(:, :), thinned(:, :)
   type(importance_sampling_result) :: object, check_is
   type(loo_result) :: result, check, collection(3)
   type(expectation_result) :: expectation, expectation_shifted
   type(comparison_result) :: comparison
   procedure(log_likelihood_callback), pointer :: callback
   integer :: i, j, m, chains(200), permutation(200)
   real(dp) :: offsets(3), d1(2), d2(2)
   character(len = 4) :: methods(3) = [character(len = 4) :: 'psis', 'sis', 'tis']
   character(len = 8) :: modes(4) = [character(len = 8) :: 'mean', 'variance', 'sd', 'quantile']
   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   inf = ieee_value(0.0_dp, ieee_positive_inf)
   eff = [0.5_dp, 1.0_dp, 0.8_dp, 1.2_dp]
   q = [0.1_dp, 0.5_dp, 0.9_dp]
   do i = 1, 200
      ll(i, 1) = -0.3_dp - 0.5_dp * sin(real(i, dp) * 0.17_dp) ** 2
      ll(i, 2) = -1.0_dp - 0.7_dp * cos(real(i, dp) * 0.113_dp) - 0.002_dp * real(i, dp)
      ll(i, 3) = -0.015_dp * real(i, dp) - 0.0004_dp * real(i, dp) ** 2
      ll(i, 4) = -2.0_dp + 0.05_dp * sin(real(i, dp) * 0.031_dp)
      lp(i) = -0.002_dp * real(i, dp)
      lg(i) = -0.003_dp * real(i, dp) + 0.05_dp * sin(real(i, dp) * 0.31_dp)
      x(i, :) = sin(real(i, dp) * 0.17_dp) + [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
      chains(i) = (i - 1) / 50 + 1
      permutation(i) = mod(i - 1, 4) * 50 + (i - 1) / 4 + 1
   end do
   a = reshape(ll, shape(a))
   do m = 1, 3
      object = importance_sampling(-ll, trim(methods(m)), eff)
      call require(object%status == 0, 'importance sampling status')
      w = importance_weights(object)
      call reference(trim(methods(m)) // '_lw', reshape(w, [size(w)]))
      call reference(trim(methods(m)) // '_k', object%pareto_k)
      call reference(trim(methods(m)) // '_ess', object%n_eff)
      call check_close(sum(exp(w), dim = 1), [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], 'normalized weights')
      result = loo(ll, eff, trim(methods(m)))
      call require(result%status == 0, 'loo status')
      call reference(trim(methods(m)) // '_pointwise', reshape(result%pointwise, [size(result%pointwise)]))
      call reference(trim(methods(m)) // '_estimates', reshape(result%estimates, [size(result%estimates)]))
      check = loo(a, eff, trim(methods(m)))
      call check_close(reshape(check%pointwise, [20]), reshape(result%pointwise, [20]), 'array layout')
      check = loo(ll - 2000.0_dp, eff, trim(methods(m)))
      call check_close(check%pointwise(:, 2), result%pointwise(:, 2), 'MCSE underflow shift', 1.0e-10_dp)
      call check_close(check%pointwise(:, 3), result%pointwise(:, 3), 'effective parameters shift', 1.0e-10_dp)
      call check_close(check%pointwise(:, 1) + 2000.0_dp, result%pointwise(:, 1), 'ELPD shift', 1.0e-10_dp)
      check = loo(ll + 1000.0_dp, eff, trim(methods(m)))
      call check_close(check%pointwise(:, 2), result%pointwise(:, 2), 'MCSE overflow shift', 1.0e-10_dp)
   end do
   object = psis(-ll, eff)
   check_is = psis(reshape(-ll, shape(a)), eff)
   call check_close(object%pareto_k, check_is%pareto_k, 'psis array')
   check_is = tis(reshape(-ll, shape(a)), eff)
   call require(check_is%status == 0, 'tis array')
   check_is = sis(reshape(-ll, shape(a)), eff)
   call require(check_is%status == 0, 'sis array')
   do j = 1, 4
      expectation = E_loo(x, object, trim(modes(j)), q, -ll)
      call require(expectation%status == 0, 'expectation status')
      call reference('expect_' // trim(modes(j)), reshape(expectation%value, [size(expectation%value)]))
      call reference('expect_k_' // trim(modes(j)), expectation%pareto_k)
   end do
   expectation = E_loo(x, object, 'variance')
   expectation_shifted = E_loo(x + 1.0e9_dp, object, 'variance')
   call check_close(reshape(expectation_shifted%value, [4]), reshape(expectation%value, [4]), &
      'weighted variance offset stability', 2.0e-7_dp)
   result = waic(ll)
   call reference('waic_pointwise', reshape(result%pointwise, [12]))
   call reference('waic_estimates', reshape(result%estimates, [6]))
   check = waic(a)
   call check_close(reshape(check%pointwise, [12]), reshape(result%pointwise, [12]), 'waic array')
   result = elpd(ll)
   call reference('elpd_pointwise', reshape(result%pointwise, [8]))
   call reference('elpd_estimates', reshape(result%estimates, [4]))
   check = elpd(a)
   call check_close(reshape(check%pointwise, [8]), reshape(result%pointwise, [8]), 'elpd array')
   result = loo_approximate_posterior(ll, lp, lg)
   call reference('ap_pointwise', reshape(result%pointwise, [20]))
   call reference('ap_estimates', reshape(result%estimates, [6]))
   check = loo_approximate_posterior(a, lp, lg)
   call check_close(reshape(check%pointwise, [20]), reshape(result%pointwise, [20]), 'approximate posterior array')
   call reference('eff', relative_eff(exp(a)))
   call reference('eff', relative_eff(exp(ll), chains))
   call reference('eff', relative_eff(exp(ll(permutation, :)), chains(permutation)))
   object = psislw(-ll)
   w = importance_weights(object)
   call reference('old_lw', reshape(w, [size(w)]))
   call reference('old_k', object%pareto_k)
   result = loo(ll, eff)
   values(1) = mcse_loo(result, threshold = 10.0_dp)
   call reference('mcse', values)
   call require(ieee_is_nan(mcse_loo(result, threshold = -10.0_dp)), 'k gating MCSE')
   offsets = [0.0_dp, 0.1_dp, 0.3_dp]
   do j = 1, 3
      do i = 1, 200
         shifted(i, :) = ll(i, :) + offsets(j) * cos(real(i, dp) * 0.29_dp)
      end do
      collection(j) = loo(shifted, eff)
      c(:, j) = collection(j)%pointwise(:, 1)
   end do
   comparison = loo_compare(c)
   call reference('compare_diff', comparison%elpd_diff)
   call reference('compare_se', comparison%se_diff)
   call reference('compare_p', comparison%p_worse)
   d1 = compare(c(:, 1), c(:, 2))
   d2 = compare(c(:, 2), c(:, 1))
   call check_close(d1, [-d2(1), d2(2)], 'comparison symmetry')
   thinned = thin_draws(reshape([(real(i, dp), i = 1,10)], [10, 1]), 3)
   call reference('thin', reshape(thinned, [3]))
   call check_close(compute_point_estimate(reshape([1.0_dp,3.0_dp,5.0_dp,7.0_dp], [2,2])), &
      [2.0_dp,6.0_dp], 'point estimate')
   callback => loglik
   check = loo_function(callback, 200, 4, eff)
   call check_close(reshape(check%pointwise, [20]), reshape(result%pointwise, [20]), 'loo callback')
   check = loo_i(3, callback, 200, [eff(3)])
   call check_close(check%pointwise(1, :), result%pointwise(3, :), 'single callback')
   check = waic_function(callback, 200, 4)
   call reference('waic_pointwise', reshape(check%pointwise, [12]))
   check = loo_approximate_posterior_function(callback, 200, 4, lp, lg)
   call reference('ap_pointwise', reshape(check%pointwise, [20]))
   callback => likelihood
   call reference('eff', relative_eff_function(callback, 200, 4, chains))
   nullify(callback)
   check = loo_function(callback, 200, 4)
   call require(check%status == 1, 'null callback rejected')
   call require(all(ieee_is_nan(relative_eff_function(callback, 200, 4, chains))), 'null likelihood callback')
   ! Degenerate and forbidden values exercise nonfinite handling without invalid arithmetic.
   t = 1.0e20_dp
   object = psis(t)
   w = importance_weights(object, logarithm = .false.)
   call check_close(w(:, 1), spread(0.01_dp, 1, 100), 'huge constant log ratios normalize')
   call check_close(object%n_eff, [100.0_dp], 'huge constant ESS')
   call require(object%pareto_k(1) == inf, 'constant-tail k infinity')
   t = ieee_value(0.0_dp, ieee_negative_inf)
   t(100) = 0.0_dp
   object = psis(t)
   call require(object%status == 0, 'negative infinity is a zero ratio')
   call check_close(object%n_eff, [1.0_dp], 'single surviving ratio ESS')
   t = ieee_value(0.0_dp, ieee_negative_inf)
   object = psis(t)
   call require(object%status == 1, 'all zero ratios rejected')
   t = 0.0_dp
   t(1) = nan
   object = psis(t)
   call require(object%status == 1, 'NaN rejected')
   t(1) = inf
   object = tis(t)
   call require(object%status == 1, 'positive infinite ratio rejected')
   object = sis([0.0_dp,1.0_dp], [-1.0_dp])
   call require(object%status == 1, 'negative efficiency rejected')
   object = sis([0.0_dp,1.0_dp], [nan])
   call require(object%status == 0, 'NaN efficiency becomes one')
   shifted = ll
   shifted(1, 1) = -inf
   check = loo(shifted)
   call require(check%status == 1, 'nonfinite likelihood rejected')
   check = waic(shifted)
   call require(check%status == 1, 'nonfinite WAIC likelihood rejected')
   chains(1) = 0
   call require(all(ieee_is_nan(relative_eff(exp(ll), chains))), 'zero chain IDs rejected')
   print *, 'Core actual-R references, shift identities, callbacks, and edge cases passed.'
contains
   subroutine loglik(observation, output)
      integer, intent(in) :: observation !! One-based observation column requested by the callback test.
      real(dp), intent(out) :: output(:) !! Log likelihood vector in original posterior draw order.
      output = ll(:, observation)
   end subroutine loglik
   subroutine likelihood(observation, output)
      integer, intent(in) :: observation !! One-based observation column requested by relative-efficiency evaluation.
      real(dp), intent(out) :: output(:) !! Unlogged likelihood vector in posterior draw order.
      output = exp(ll(:, observation))
   end subroutine likelihood
   subroutine require(condition, label)
      logical, intent(in) :: condition !! Assertion that must hold for the test to pass.
      character(*), intent(in) :: label !! Description printed on assertion failure.
      if (.not. condition) then
         print *, 'FAILED: ', label
         error stop 1
      end if
   end subroutine require
   subroutine check_close(actual, expected, label, tolerance)
      real(dp), intent(in) :: actual(:) !! Numerically computed test results.
      real(dp), intent(in) :: expected(:) !! Expected values from R or an independent identity.
      character(*), intent(in) :: label !! Description of the equality being checked.
      real(dp), intent(in), optional :: tolerance !! Relative and absolute tolerance, default 2e-10.
      real(dp) :: tol
      integer :: k
      tol = 2.0e-10_dp
      if (present(tolerance)) tol = tolerance
      call require(size(actual) == size(expected), label // ' size')
      do k = 1, size(actual)
         if (ieee_is_nan(expected(k))) then
            call require(ieee_is_nan(actual(k)), label // ' NaN')
         else if (abs(actual(k) - expected(k)) > tol * max(1.0_dp, abs(expected(k))) .or. ieee_is_nan(actual(k))) then
            print *, label, k, ' actual, expected: ', actual(k), expected(k)
            error stop 2
         end if
      end do
   end subroutine check_close
   subroutine reference(name, actual)
      character(*), intent(in) :: name !! Base name of the actual-R fixture under test/reference/data.
      real(dp), intent(in) :: actual(:) !! Values to compare, flattened in Fortran/R column-major order.
      real(dp) :: expected(size(actual))
      integer :: unit, k
      open(newunit = unit, file = 'test/reference/data/' // name // '.txt', status = 'old', action = 'read')
      do k = 1, size(actual)
         read(unit, *) expected(k)
      end do
      close(unit)
      call check_close(actual, expected, name)
   end subroutine reference
end program test_core
