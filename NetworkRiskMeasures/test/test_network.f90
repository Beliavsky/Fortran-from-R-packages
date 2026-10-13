! SPDX-License-Identifier: GPL-3.0-only
program test_network
   use network_risk_measures
   use network_reference_values
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan, ieee_positive_inf
   implicit none
   real(dp) :: a(3, 3), b(3), w(3), m(3, 3), v(4, 4), c(4, 4), c2(4, 4), shocks(4, 3)
   real(dp) :: s(3), f, d(3, 3), rows(7), cols(7), estimate(7, 7), other(7, 7), error, uniforms(10000)
   real(dp) :: bad(3, 3), singleton(1, 1), sumry(5, 1), two(2, 2), line(5, 5), output(5, 5)
   real(dp) :: one_shock(5, 1), expected(5), negative(3), inf
   type(contagion_result) :: fit, alternative
   integer :: status, i, j, steps, draws, tests
   tests = 0
   a = reshape([0.0_dp, 10.0_dp, 3.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 3.0_dp, 0.0_dp], [3, 3])
   b = [2.0_dp, 5.0_dp, 2.0_dp]
   w = [10.0_dp, 100.0_dp, 30.0_dp]
   call risk_matrix(a, m, status, b)
   call good(status, 'risk assets')
   call near(reshape(m, [9]), ref_risk_assets, 1.0e-14_dp, 'R risk assets')
   call risk_matrix(a, m, status, b, exposure_type='liabilities')
   call near(reshape(m, [9]), ref_risk_liabilities, 1.0e-14_dp, 'R risk liabilities')
   call vulnerability_matrix(a, m, status, b)
   call near(reshape(m, [9]), ref_vulnerability, 1.0e-14_dp, 'R vulnerability')
   call impact_matrix(transpose(a), m, status, b, exposure_type='liabilities')
   call near(reshape(m, [9]), ref_risk_assets, 1.0e-14_dp, 'transpose convention')
   call criticality(a, s, status, b)
   call near(s, ref_criticality, 1.0e-14_dp, 'R criticality')
   v = reshape([0.0_dp, 0.2_dp, 0.1_dp, 0.0_dp, 0.7_dp, 0.0_dp, 0.4_dp, 0.0_dp, &
      0.0_dp, 0.6_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [4, 4])
   call communicability_matrix(v, c, status)
   call good(status, 'full exponential')
   call near(reshape(c, [16]), ref_comm_full, 2.0e-12_dp, 'R full communicability')
   call communicability_matrix(v, c, status, terms=3)
   call near(reshape(c, [16]), ref_comm_three, 2.0e-14_dp, 'R truncated communicability')
   call communicability_matrix(v, c, status, terms=0)
   c2 = 0.0_dp
   do i = 1, 4
      c2(i, i) = 1.0_dp
   end do
   call near(reshape(c, [16]), reshape(c2, [16]), 0.0_dp, 'zero degree identity')
   call communicability_matrix(v, c, status, terms=200)
   call near(reshape(c, [16]), ref_comm_full, 2.0e-12_dp, 'high-degree factorial-safe recurrence')
   call impact_susceptibility(a, s, status, b)
   call near(s, ref_susceptibility, 2.0e-12_dp, 'R susceptibility')
   call impact_susceptibility(a, s, status, b, weights=w, binary=.false.)
   call near(s, ref_weighted_susceptibility, 2.0e-11_dp, 'R weighted susceptibility')
   call impact_fluidity(a, f, status, b)
   call near([f], ref_fluidity, 2.0e-12_dp, 'R fluidity')
   call impact_fluidity(a, f, status, b, weights=w, binary=.false., terms=7)
   call near([f], ref_weighted_fluidity, 2.0e-11_dp, 'R weighted fluidity')
   call impact_diffusion(a, d, status, b)
   call near(reshape(d, [9]), ref_diffusion, 2.0e-12_dp, 'R diffusion')
   call impact_diffusion(a, d, status, b, weights=w, binary=.false.)
   call near(reshape(d, [9]), ref_weighted_diffusion, 2.0e-11_dp, 'R weighted diffusion')
   shocks(:, 1) = [0.4_dp, 0.0_dp, 0.0_dp, 0.0_dp]
   shocks(:, 2) = [0.0_dp, 0.6_dp, 0.0_dp, 0.1_dp]
   shocks(:, 3) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
   call contagion(v, fit, status, buffer=[2.0_dp, 3.0_dp, 5.0_dp, 7.0_dp], shocks=shocks, &
      weights=[1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], exposure_type='impact', abs_tol=1.0e-12_dp, max_it=1000, keep_history=.true.)
   call good(status, 'DebtRank convergence')
   call near(reshape(fit%final, [12]), ref_debtrank_final, 2.0e-12_dp, 'R DebtRank final')
   call near(reshape(fit%summary, [15]), ref_debtrank_summary, 2.0e-11_dp, 'R DebtRank summary')
   call near(real(fit%iterations, dp), ref_debtrank_iterations, 0.0_dp, 'R DebtRank iterations')
   do i = 1, 3
      call near(fit%history(:, fit%iterations(i), i), fit%final(:, i), 0.0_dp, 'last history state')
      call check(all(ieee_is_nan(fit%history(:, fit%iterations(i) + 1:, i))), 'unused history NaN')
   end do
   call contagion(v, fit, status, buffer=[2.0_dp, 3.0_dp, 5.0_dp, 7.0_dp], shocks=shocks, &
      weights=[1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], exposure_type='impact', abs_tol=1.0e-12_dp, &
      max_it=1000, single_hit=.true.)
   call near(reshape(fit%final, [12]), ref_single_final, 2.0e-12_dp, 'R single-hit final')
   call near(reshape(fit%summary, [15]), ref_single_summary, 2.0e-11_dp, 'R single-hit summary')
   call near(real(fit%iterations, dp), ref_single_iterations, 0.0_dp, 'R single-hit iterations')
   call check(.not. allocated(fit%history), 'history omitted by default')
   call contagion(v, fit, status, buffer=[2.0_dp, 3.0_dp, 5.0_dp, 7.0_dp], shocks=shocks, &
      weights=[1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], exposure_type='impact', abs_tol=1.0e-12_dp, &
      max_it=1000, method='threshold', alpha=0.5_dp)
   call near(reshape(fit%final, [12]), ref_threshold_final, 2.0e-12_dp, 'R threshold final')
   call near(reshape(fit%summary, [15]), ref_threshold_summary, 2.0e-11_dp, 'R threshold summary')
   call near(real(fit%iterations, dp), ref_threshold_iterations, 0.0_dp, 'R threshold iterations')
   rows = [7.0_dp, 5.0_dp, 3.0_dp, 1.0_dp, 3.0_dp, 0.0_dp, 1.0_dp]
   cols = [4.0_dp, 5.0_dp, 5.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, 4.0_dp]
   call max_ent(rows, cols, estimate, status, abs_tol=1.0e-10_dp, iterations=steps, residual=error)
   call good(status, 'maximum entropy converged')
   call near(reshape(estimate, [49]), ref_max_ent, 1.0e-12_dp, 'R maximum entropy')
   call check(error < 1.0e-10_dp .and. steps > 1, 'maximum entropy diagnostics')
   call matrix_estimation(rows, cols, other, status, abs_tol=1.0e-10_dp)
   call near(reshape(other, [49]), ref_matrix_estimation, 1.0e-12_dp, 'R matrix estimation dispatcher')
   do i = 1, size(uniforms)
      uniforms(i) = real(modulo(i * 37, 997), dp) / 997.0_dp
   end do
   call min_dens(rows, cols, estimate, status, cost=0.3_dp, lambda=0.7_dp, k=20, remove_prob=0.08_dp, &
      abs_tol=1.0e-10_dp, max_it=10000, uniforms=uniforms, iterations=steps, residual=error, draws_used=draws)
   call good(status, 'controlled minimum density converged')
   call near(reshape(estimate, [49]), ref_min_dens_controlled, 2.0e-12_dp, 'R controlled-draw minimum density')
   call near([real(draws, dp)], ref_min_dens_draws, 0.0_dp, 'R controlled uniform count')
   call check(error < 2.0e-10_dp, 'minimum density residual')
   call min_dens(rows, cols, estimate, status, seed=192, max_it=100000)
   call good(status, 'seeded minimum density')
   call min_dens(rows, cols, other, status, seed=192, max_it=100000)
   call near(reshape(estimate, [49]), reshape(other, [49]), 0.0_dp, 'local RNG repeatability')
   call matrix_estimation(rows, cols, other, status, method='md', seed=192)
   call near(reshape(estimate, [49]), reshape(other, [49]), 0.0_dp, 'minimum-density dispatcher')
   call near(sum(estimate, dim=2), rows, 1.0e-3_dp, 'seeded row margins')
   call near(sum(estimate, dim=1), cols, 1.0e-3_dp, 'seeded column margins')
   do i = 1, 7
      call check(estimate(i, i) == 0.0_dp, 'zero reconstruction diagonal')
   end do
   do i = 1, 10
      call min_dens(rows, cols, estimate, status, seed=i, abs_tol=1.0e-10_dp, residual=error)
      call good(status, 'minimum-density multi-seed stress convergence')
      call check(error < 2.0e-10_dp .and. all(estimate >= 0.0_dp), 'minimum-density stress margins and sign')
   end do
   ! Independent exact directed-chain formulas and orientation checks.
   line = 0.0_dp
   do i = 1, 4
      line(i, i + 1) = 0.5_dp
   end do
   call communicability_matrix(line, output, status)
   do j = 1, 5
      do i = 1, 5
         f = 0.0_dp
         if (j >= i) then
            f = 1.0_dp
            do steps = 1, j - i
               f = f * 0.5_dp / real(steps, dp)
            end do
         end if
         call near([output(i, j)], [f], 2.0e-14_dp, 'analytic directed-chain exponential')
      end do
   end do
   one_shock = 0.0_dp
   one_shock(1, 1) = 0.8_dp
   expected = [0.8_dp, 0.4_dp, 0.2_dp, 0.1_dp, 0.05_dp]
   call contagion(line, fit, status, shocks=one_shock, exposure_type='impact', abs_tol=1.0e-14_dp)
   call near(fit%final(:, 1), expected, 1.0e-14_dp, 'analytic DebtRank chain')
   call contagion(line, fit, status, shocks=one_shock, exposure_type='impact', method='threshold')
   call near(fit%final(:, 1), one_shock(:, 1), 0.0_dp, 'sub-threshold shock does not propagate')
   line(1, 2) = 2.0_dp
   one_shock(1, 1) = 1.0_dp
   call contagion(line, fit, status, shocks=one_shock, exposure_type='impact')
   call check(fit%final(2, 1) > 1.0_dp, 'uncapped final stress preserved')
   call contagion(line, fit, status, shocks=one_shock, exposure_type='impact', max_it=1)
   call check(status == nr_iterations .and. .not. fit%converged(1), 'contagion iteration exhaustion')
   call near(fit%final(:, 1), [1.0_dp, 2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], 0.0_dp, 'last nonconverged iterate')
   ! Empty-edge, singleton, shape and invalid-domain cases.
   bad = 0.0_dp
   call impact_susceptibility(bad, s, status, exposure_type='impact')
   call near(s, [0.0_dp, 0.0_dp, 0.0_dp], 0.0_dp, 'disconnected susceptibility')
   call impact_diffusion(bad, d, status, exposure_type='impact')
   call check(all(d == 0.0_dp), 'disconnected diffusion')
   call contagion(bad, fit, status, exposure_type='impact')
   call near(reshape(fit%final, [9]), reshape(fit%initial, [9]), 0.0_dp, 'disconnected shocks')
   call check(all(fit%iterations == 1), 'zero network stops on first step')
   singleton = 0.0_dp
   c(1, 1) = 0.0_dp
   call communicability_matrix(c(:1, :1), singleton, status)
   call check(status == nr_ok .and. singleton(1, 1) == 1.0_dp, 'singleton zero exponential')
   c(1, 1) = -2.0_dp
   call communicability_matrix(c(:1, :1), singleton, status)
   call near([singleton(1, 1)], [exp(-2.0_dp)], 1.0e-14_dp, 'signed scalar exponential')
   call max_ent([0.0_dp], [0.0_dp], singleton, status)
   call good(status, 'zero singleton maximum entropy')
   call check(singleton(1, 1) == 0.0_dp, 'zero singleton reconstruction')
   call min_dens([0.0_dp], [0.0_dp], singleton, status)
   call good(status, 'zero singleton minimum density')
   call max_ent([1.0_dp], [1.0_dp], singleton, status)
   call check(status == nr_infeasible, 'positive singleton infeasible')
   call max_ent(rows, cols, estimate, status, max_it=1, abs_tol=1.0e-12_dp)
   call check(status == nr_iterations, 'maximum entropy iteration exhaustion')
   call max_ent([-1.0_dp, 1.0_dp], [0.0_dp, 0.0_dp], two, status)
   call check(status == nr_domain, 'negative marginal rejected')
   call max_ent([2.0_dp, 1.0_dp], [1.0_dp, 1.0_dp], two, status)
   call check(status == nr_infeasible, 'unequal totals rejected')
   call max_ent([2.0_dp, 1.0_dp], [2.0_dp, 1.0_dp], two, status)
   call check(status == nr_infeasible, 'diagonal-only margin conflict rejected')
   call max_ent([2.0_dp, 3.0_dp], [3.0_dp, 2.0_dp], two, status, abs_tol=1.0e-12_dp)
   call good(status, 'exact two-vertex maximum entropy')
   call near(reshape(two, [4]), [0.0_dp, 3.0_dp, 2.0_dp, 0.0_dp], 0.0_dp, 'analytic two-vertex margins')
   call min_dens([2.0_dp, 3.0_dp], [3.0_dp, 2.0_dp], two, status, cost=0.0_dp, remove_prob=1.0_dp, &
      max_it=2, uniforms=[0.5_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.0_dp], draws_used=draws)
   call check(status == nr_iterations .and. all(two == 0.0_dp) .and. draws == 5, 'singleton occupied removal')
   call min_dens(rows, cols, estimate, status, uniforms=uniforms(:1))
   call check(status == nr_random_exhausted, 'uniform stream exhaustion')
   call min_dens(rows, cols, estimate, status, uniforms=[1.0_dp])
   call check(status == nr_domain, 'out-of-range uniform rejected')
   call min_dens(rows, cols, estimate, status, lambda=-0.1_dp)
   call check(status == nr_domain, 'negative lambda rejected')
   call min_dens(rows, cols, estimate, status, max_it=1, cost=1000.0_dp, uniforms=[0.9_dp, 0.5_dp, 0.5_dp])
   call check(status == nr_iterations .and. all(estimate == 0.0_dp), 'rejected Metropolis proposal')
   call min_dens(rows, cols, estimate, status, max_it=1, cost=0.0_dp, uniforms=[0.9_dp, 0.5_dp, 0.5_dp])
   call check(status == nr_iterations .and. sum(estimate) > 0.0_dp, 'accepted Metropolis proposal')
   call risk_matrix(a, m, status, buffer=[0.0_dp, 1.0_dp, 1.0_dp])
   call check(status == nr_domain .and. all(ieee_is_nan(m)), 'zero buffer rejected')
   negative = [-1.0_dp, 1.0_dp, 1.0_dp]
   call risk_matrix(a, m, status, buffer=negative)
   call check(status == nr_domain, 'negative buffer rejected')
   call risk_matrix(a, m, status, buffer=b(:2))
   call check(status == nr_shape, 'short buffer rejected')
   call risk_matrix(a(:, :2), m, status, buffer=b)
   call check(status == nr_shape, 'nonsquare network rejected')
   bad = a
   bad(1, 1) = -1.0_dp
   call risk_matrix(bad, m, status, b)
   call check(status == nr_domain, 'negative network rejected')
   bad(1, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
   call risk_matrix(bad, m, status, b)
   call check(status == nr_domain, 'NaN network rejected')
   inf = ieee_value(0.0_dp, ieee_positive_inf)
   bad(1, 1) = inf
   call communicability_matrix(bad, m, status)
   call check(status == nr_domain, 'infinite exponential input rejected')
   call communicability_matrix(a, m, status, terms=-2)
   call check(status == nr_domain, 'negative Taylor degree rejected')
   call contagion(a, fit, status, b, weights=[0.0_dp, 0.0_dp, 0.0_dp])
   call check(status == nr_domain, 'zero weights rejected')
   call contagion(a, fit, status, b, weights=[huge(1.0_dp), huge(1.0_dp), huge(1.0_dp)])
   call good(status, 'stable weight normalization')
   call contagion(a, alternative, status, b)
   call near(reshape(fit%summary, [15]), reshape(alternative%summary, [15]), 1.0e-13_dp, 'weight scale invariance')
   call contagion(a, fit, status, b, method='missing')
   call check(status == nr_domain, 'unknown propagation rejected')
   call contagion(a, fit, status, b, alpha=0.0_dp)
   call check(status == nr_domain, 'invalid propagation trigger rejected')
   call summarize_contagion(reshape([0.0_dp], [1, 1]), reshape([2.0_dp], [1, 1]), &
      [1.0_dp], [5.0_dp], sumry, status, cap=0.5_dp)
   call near(sumry(:, 1), [0.0_dp, 0.5_dp, 0.0_dp, 10.0_dp, 1.0_dp], 0.0_dp, 'summary cap leaves losses raw')
   ! Deterministic dense stress/property checks over twenty varying networks.
   do steps = 1, 20
      do j = 1, 4
         do i = 1, 4
            v(i, j) = real(modulo(i * 7 + j * 11 + steps * 13, 31), dp) / 100.0_dp
         end do
         v(j, j) = 0.0_dp
      end do
      call communicability_matrix(v, c, status)
      call communicability_matrix(transpose(v), c2, status)
      call near(reshape(c, [16]), reshape(transpose(c2), [16]), 2.0e-12_dp, 'exponential transpose property')
      call contagion(v, fit, status, exposure_type='impact', abs_tol=1.0e-12_dp, max_it=1000)
      call good(status, 'dense stress convergence')
      call check(all(fit%final >= fit%initial), 'nonnegative propagation monotonicity')
      call contagion(transpose(v), alternative, status, exposure_type='vulnerability', abs_tol=1.0e-12_dp, max_it=1000)
      call near(reshape(fit%final, [16]), reshape(alternative%final, [16]), 1.0e-12_dp, 'contagion orientation property')
   end do
   print '(a,i0)', 'NetworkRiskMeasures checks passed: ', tests
contains
   subroutine check(condition, label)
      logical, intent(in) :: condition !! Assertion predicate that must be true for this deterministic test.
      character(*), intent(in) :: label !! Description printed before stopping if the assertion fails.
      tests = tests + 1
      if (.not. condition) then
         print '(a)', 'FAILED: ' // label
         error stop 1
      end if
   end subroutine check
   subroutine good(code, label)
      integer, intent(in) :: code !! Public routine status expected to equal nr_ok.
      character(*), intent(in) :: label !! Human-readable operation whose status is being checked.
      if (code /= nr_ok) print '(a,i0)', 'Unexpected status: ', code
      call check(code == nr_ok, label)
   end subroutine good
   subroutine near(actual, wanted, tolerance, label)
      real(dp), intent(in) :: actual(:) !! Computed numerical values to compare elementwise.
      real(dp), intent(in) :: wanted(:) !! Reference values of the same length.
      real(dp), intent(in) :: tolerance !! Absolute nonnegative permitted maximum discrepancy.
      character(*), intent(in) :: label !! Explanation of the reference or independent property being checked.
      call check(size(actual) == size(wanted), label // ' shape')
      if (any(ieee_is_nan(actual))) then
         call check(.false., label // ' NaN')
      end if
      if (maxval(abs(actual - wanted)) > tolerance) then
         print '(a,es24.16)', 'Maximum discrepancy: ', maxval(abs(actual - wanted))
      end if
      call check(maxval(abs(actual - wanted)) <= tolerance, label)
   end subroutine near
end program test_network
