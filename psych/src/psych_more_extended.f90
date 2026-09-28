module psych_more_extended
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use r_kinds, only : dp
   use r_linalg, only : symmetric_eigen
   use psych_correlations, only : pseudoinverse, squared_multiple_correlations
   use psych_data_extended, only : mixed_correlation_matrix, stats_by_numeric
   use psych_factor, only : promax_rotate, procrustes_rotate
   use psych_factor_extended, only : minimum_residual_factor, eigenvalue_factor_count, factor_extension_loadings
   use psych_iclust_extended, only : iclust_hierarchical
   use psych_omega_extended, only : omega_from_correlation
   use psych_scoring_extended, only : score_overlap_summary
   use psych_simulation, only : simulate_correlated_normal, simulate_npl_items, simulate_npn_items
   use psych_computational_extended, only : irt_factor_analysis
   use psych_irt, only : two_parameter_person_scores
   use psych_validity_extended, only : yule_correlation_matrix, phi_to_tetrachoric_values
   use psych_simulation_algorithms, only : simulate_minor_factors
   use psych_stats, only : correlation_matrix
   implicit none
   private
   public :: cta_dynamics, cta15_dynamics, stats_by_bootstrap, score_by_correlations
   public :: mediate_from_covariance, factor_analysis_sapa, factor_analysis_pooled
   public :: parallel_analysis_poly, simulate_omega_study
   public :: bass_ackward_factor, factor_analysis_by_groups, simulate_multi_subject
   public :: best_scales_from_correlation, exploratory_sem, correlation_fiml
   public :: test_irt_benchmark, phi_demo_simulation, glb_algebraic_admm

contains

   subroutine cta15_dynamics(cues, n_steps, tendency_time, action_time, time_fraction, mean_tendency, mean_action, &
      stimulus, inhibition, consumption, initial_tendency, initial_action, info)
      real(dp), intent(in) :: cues(:) !! Fixed cue-strength vector that excites action tendencies.
      integer, intent(in) :: n_steps !! Number of discrete 0.05-step iterations to simulate; must be positive.
      real(dp), allocatable, intent(out) :: tendency_time(:, :) !! Step-by-action tendency trajectory.
      real(dp), allocatable, intent(out) :: action_time(:, :) !! Step-by-action nonnegative action-strength trajectory.
      real(dp), allocatable, intent(out) :: time_fraction(:) !! Fraction of iterations won by each action.
      real(dp), allocatable, intent(out) :: mean_tendency(:) !! Mean tendency over all simulated iterations.
      real(dp), allocatable, intent(out) :: mean_action(:) !! Mean action strength over all simulated iterations.
      real(dp), intent(in), optional :: stimulus(:, :) !! Cue-to-tendency stimulation matrix; identity by default.
      real(dp), intent(in), optional :: inhibition(:, :) !! Action inhibition matrix; ones with 0.05 diagonal by default.
      real(dp), intent(in), optional :: consumption(:, :) !! Action-to-tendency consumption matrix; 0.05 identity by default.
      real(dp), intent(in), optional :: initial_tendency(:) !! Initial tendency vector; zeros by default.
      real(dp), intent(in), optional :: initial_action(:) !! Initial action-strength vector; cue strengths by default.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid sizes or iteration count.
      real(dp), allocatable :: stim(:, :), inhibit(:, :), consume(:, :), tendency(:), action(:)
      integer, allocatable :: counts(:)
      integer :: i, j, n, winner, local_info
      real(dp), parameter :: step = 0.05_dp

      n = size(cues)
      local_info = 0
      if (n < 1 .or. n_steps < 1) local_info = 1
      if (present(stimulus)) then
         if (size(stimulus, 1) /= n .or. size(stimulus, 2) /= n) local_info = 2
      end if
      if (present(inhibition)) then
         if (size(inhibition, 1) /= n .or. size(inhibition, 2) /= n) local_info = 3
      end if
      if (present(consumption)) then
         if (size(consumption, 1) /= n .or. size(consumption, 2) /= n) local_info = 4
      end if
      if (present(initial_tendency)) then
         if (size(initial_tendency) /= n) local_info = 5
      end if
      if (present(initial_action)) then
         if (size(initial_action) /= n) local_info = 6
      end if
      if (local_info /= 0) then
         allocate(tendency_time(0, 0), action_time(0, 0), time_fraction(0), mean_tendency(0), mean_action(0))
         if (present(info)) info = local_info
         return
      end if

      allocate(stim(n, n), inhibit(n, n), consume(n, n), tendency(n), action(n), counts(n))
      stim = 0.0_dp
      inhibit = 1.0_dp
      consume = 0.0_dp
      do i = 1, n
         stim(i, i) = 1.0_dp
         inhibit(i, i) = 0.05_dp
         consume(i, i) = 0.05_dp
      end do
      if (present(stimulus)) stim = stimulus
      if (present(inhibition)) inhibit = inhibition
      if (present(consumption)) consume = consumption
      tendency = 0.0_dp
      action = cues
      if (present(initial_tendency)) tendency = initial_tendency
      if (present(initial_action)) action = initial_action

      allocate(tendency_time(n_steps, n), action_time(n_steps, n), time_fraction(n), mean_tendency(n), mean_action(n))
      counts = 0
      mean_tendency = 0.0_dp
      mean_action = 0.0_dp
      do i = 1, n_steps
         tendency = tendency + step * matmul(cues, stim) - step * matmul(action, consume)
         action = action + step * tendency - step * matmul(action, inhibit)
         do j = 1, n
            if (action(j) < 0.0_dp) action(j) = 0.0_dp
         end do
         tendency_time(i, :) = tendency
         action_time(i, :) = action
         mean_tendency = mean_tendency + tendency
         mean_action = mean_action + action
         winner = maxloc(action, dim=1)
         counts(winner) = counts(winner) + 1
      end do
      mean_tendency = mean_tendency / real(n_steps, dp)
      mean_action = mean_action / real(n_steps, dp)
      time_fraction = real(counts, dp) / real(n_steps, dp)
      if (present(info)) info = 0
   end subroutine cta15_dynamics

   subroutine cta_dynamics(cues, n_steps, tendency_time, action_time, time_fraction, final_expectation, &
      inhibition, expectation, consumption, tendency_strength, reward, initial_tendency, initial_action, learn, info)
      real(dp), intent(in) :: cues(:) !! Fixed cue-strength vector that excites action tendencies.
      integer, intent(in) :: n_steps !! Number of discrete 0.05-step iterations to simulate; must be positive.
      real(dp), allocatable, intent(out) :: tendency_time(:, :) !! Step-by-action tendency trajectory.
      real(dp), allocatable, intent(out) :: action_time(:, :) !! Step-by-action nonnegative action-strength trajectory.
      real(dp), allocatable, intent(out) :: time_fraction(:) !! Fraction of iterations won by each action.
      real(dp), allocatable, intent(out) :: final_expectation(:, :) !! Final cue-to-outcome expectation matrix.
      real(dp), intent(in), optional :: inhibition(:, :) !! Action inhibition matrix; ones with 0.05 diagonal by default.
      real(dp), intent(in), optional :: expectation(:, :) !! Initial cue-to-outcome expectation matrix; identity by default.
      real(dp), intent(in), optional :: consumption(:, :) !! Action-to-tendency consumption matrix; 0.03 identity by default.
      real(dp), intent(in), optional :: tendency_strength(:, :) !! Tendency-to-action strength matrix; identity by default.
      real(dp), intent(in), optional :: reward(:) !! Reward on each winning-action diagonal; default last action 0.01.
      real(dp), intent(in), optional :: initial_tendency(:) !! Initial tendency vector; zeros by default.
      real(dp), intent(in), optional :: initial_action(:) !! Initial action-strength vector; cue strengths by default.
      logical, intent(in), optional :: learn !! Update expectations on winner transitions when true; default true.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid dimensions or iteration count.
      real(dp), allocatable :: inhibit(:, :), expect(:, :), consume(:, :), tstrength(:, :), rew(:)
      real(dp), allocatable :: tendency(:), action(:), diag_update(:)
      integer, allocatable :: counts(:)
      real(dp) :: trace_expect
      integer :: i, j, n, winner, old_winner, local_info
      logical :: do_learn
      real(dp), parameter :: step = 0.05_dp

      n = size(cues)
      local_info = 0
      if (n < 1 .or. n_steps < 1) local_info = 1
      if (present(inhibition)) then
         if (size(inhibition, 1) /= n .or. size(inhibition, 2) /= n) local_info = 2
      end if
      if (present(expectation)) then
         if (size(expectation, 1) /= n .or. size(expectation, 2) /= n) local_info = 3
      end if
      if (present(consumption)) then
         if (size(consumption, 1) /= n .or. size(consumption, 2) /= n) local_info = 4
      end if
      if (present(tendency_strength)) then
         if (size(tendency_strength, 1) /= n .or. size(tendency_strength, 2) /= n) local_info = 5
      end if
      if (present(reward)) then
         if (size(reward) /= n) local_info = 6
      end if
      if (present(initial_tendency)) then
         if (size(initial_tendency) /= n) local_info = 7
      end if
      if (present(initial_action)) then
         if (size(initial_action) /= n) local_info = 8
      end if
      if (local_info /= 0) then
         allocate(tendency_time(0, 0), action_time(0, 0), time_fraction(0), final_expectation(0, 0))
         if (present(info)) info = local_info
         return
      end if

      allocate(inhibit(n, n), expect(n, n), consume(n, n), tstrength(n, n), rew(n))
      allocate(tendency(n), action(n), diag_update(n), counts(n))
      inhibit = 1.0_dp
      expect = 0.0_dp
      consume = 0.0_dp
      tstrength = 0.0_dp
      rew = 0.0_dp
      do i = 1, n
         inhibit(i, i) = 0.05_dp
         expect(i, i) = 1.0_dp
         consume(i, i) = 0.03_dp
         tstrength(i, i) = 1.0_dp
      end do
      rew(n) = 0.01_dp
      if (present(inhibition)) inhibit = inhibition
      if (present(expectation)) expect = expectation
      if (present(consumption)) consume = consumption
      if (present(tendency_strength)) tstrength = tendency_strength
      if (present(reward)) rew = reward
      tendency = 0.0_dp
      action = cues
      if (present(initial_tendency)) tendency = initial_tendency
      if (present(initial_action)) action = initial_action
      do_learn = .true.
      if (present(learn)) do_learn = learn
      old_winner = maxloc(action, dim=1)
      counts = 0

      allocate(tendency_time(n_steps, n), action_time(n_steps, n), time_fraction(n))
      do i = 1, n_steps
         tendency = tendency + step * matmul(cues, expect) - step * matmul(action, consume)
         action = action + step * matmul(tendency, tstrength) - step * matmul(action, inhibit)
         do j = 1, n
            if (action(j) < 0.0_dp) action(j) = 0.0_dp
         end do
         tendency_time(i, :) = tendency
         action_time(i, :) = action
         winner = maxloc(action, dim=1)
         counts(winner) = counts(winner) + 1
         if (do_learn .and. winner /= old_winner) then
            diag_update = action * rew
            do j = 1, n
               expect(j, j) = expect(j, j) + diag_update(j)
            end do
            trace_expect = 0.0_dp
            do j = 1, n
               trace_expect = trace_expect + expect(j, j)
            end do
            if (abs(trace_expect) > tiny(1.0_dp)) expect = expect / trace_expect
            old_winner = winner
         end if
      end do
      time_fraction = real(counts, dp) / real(n_steps, dp)
      allocate(final_expectation(n, n))
      final_expectation = expect
      if (present(info)) info = 0
   end subroutine cta_dynamics

   subroutine stats_by_bootstrap(x, groups, n_trials, icc1_trials, icc2_trials, f_trials, replace, info)
      real(dp), intent(in) :: x(:, :) !! Observation-by-variable numeric data matrix.
      integer, intent(in) :: groups(:) !! Integer group label paired with each observation.
      integer, intent(in) :: n_trials !! Number of label-resampling bootstrap/permutation trials.
      real(dp), allocatable, intent(out) :: icc1_trials(:, :) !! Trial-by-variable one-way single-measure ICC estimates.
      real(dp), allocatable, intent(out) :: icc2_trials(:, :) !! Trial-by-variable group-mean reliability estimates.
      real(dp), allocatable, intent(out) :: f_trials(:, :) !! Trial-by-variable between/within mean-square ratios.
      logical, intent(in), optional :: replace !! Sample group labels with replacement; default true.
      integer, intent(out), optional :: info !! Zero on success; nonzero for dimension or grouped-statistics failures.
      integer, allocatable :: sampled_groups(:), group_values(:), counts(:, :), permutation(:)
      real(dp), allocatable :: means(:, :), sds(:, :), grand(:), msb(:), msw(:), f_values(:), icc1(:), icc2(:)
      integer :: i, j, n, local_info, stat_info, pick, temp
      logical :: do_replace
      real(dp) :: u

      n = size(groups)
      local_info = 0
      if (size(x, 1) /= n .or. n < 2 .or. n_trials < 1) local_info = 1
      if (local_info /= 0) then
         allocate(icc1_trials(0, 0), icc2_trials(0, 0), f_trials(0, 0))
         if (present(info)) info = local_info
         return
      end if
      do_replace = .true.
      if (present(replace)) do_replace = replace
      allocate(sampled_groups(n), permutation(n))
      allocate(icc1_trials(n_trials, size(x, 2)), icc2_trials(n_trials, size(x, 2)), f_trials(n_trials, size(x, 2)))
      icc1_trials = 0.0_dp
      icc2_trials = 0.0_dp
      f_trials = 0.0_dp

      do i = 1, n_trials
         if (do_replace) then
            do j = 1, n
               call random_number(u)
               pick = min(n, int(u * real(n, dp)) + 1)
               sampled_groups(j) = groups(pick)
            end do
         else
            do j = 1, n
               permutation(j) = j
            end do
            do j = n, 2, -1
               call random_number(u)
               pick = min(j, int(u * real(j, dp)) + 1)
               temp = permutation(j)
               permutation(j) = permutation(pick)
               permutation(pick) = temp
            end do
            do j = 1, n
               sampled_groups(j) = groups(permutation(j))
            end do
         end if
         call stats_by_numeric(x, sampled_groups, group_values, means, sds, counts, grand, msb, msw, &
            f_values, icc1, icc2, stat_info)
         if (stat_info /= 0) local_info = max(local_info, 10 + stat_info)
         if (size(icc1) == size(x, 2)) icc1_trials(i, :) = icc1
         if (size(icc2) == size(x, 2)) icc2_trials(i, :) = icc2
         if (size(f_values) == size(x, 2)) f_trials(i, :) = f_values
      end do
      if (present(info)) info = local_info
   end subroutine stats_by_bootstrap

   subroutine score_by_correlations(keys, correlations, adjusted_correlations, alpha, average_r, g6, scale_variance, &
      pair_counts, min_n, info)
      real(dp), intent(in) :: keys(:, :) !! Item-by-scale signed key matrix shared across groups.
      real(dp), intent(in) :: correlations(:, :, :) !! Item-by-item-by-group correlation matrices from grouped data.
      real(dp), allocatable, intent(out) :: adjusted_correlations(:, :, :) !! Overlap-adjusted scale correlations by group.
      real(dp), allocatable, intent(out) :: alpha(:, :) !! Scale-by-group coefficient alpha values.
      real(dp), allocatable, intent(out) :: average_r(:, :) !! Scale-by-group average inter-item correlations.
      real(dp), allocatable, intent(out) :: g6(:, :) !! Scale-by-group Guttman lambda-six estimates.
      real(dp), allocatable, intent(out) :: scale_variance(:, :) !! Scale-by-group raw keyed variances.
      integer, intent(in), optional :: pair_counts(:, :, :) !! Pairwise observation counts corresponding to correlations.
      integer, intent(in), optional :: min_n !! Minimum mean pairwise count required to score a group; default three.
      integer, allocatable, intent(out), optional :: info(:) !! Per-group zero success or nonzero skip/computation status.
      real(dp), allocatable :: adj(:, :), a(:), av(:), g(:), item_cor(:, :), raw_cov(:, :)
      integer, allocatable :: local_status(:)
      integer :: grp, ng, ns, local_info, minimum_n
      real(dp) :: mean_count

      ng = size(correlations, 3)
      ns = size(keys, 2)
      minimum_n = 3
      if (present(min_n)) minimum_n = max(1, min_n)
      allocate(adjusted_correlations(ns, ns, ng), alpha(ns, ng), average_r(ns, ng), g6(ns, ng))
      allocate(scale_variance(ns, ng), local_status(ng))
      adjusted_correlations = 0.0_dp
      alpha = 0.0_dp
      average_r = 0.0_dp
      g6 = 0.0_dp
      scale_variance = 0.0_dp
      local_status = 0
      if (size(correlations, 1) /= size(keys, 1) .or. size(correlations, 2) /= size(keys, 1)) then
         local_status = 1
         if (present(info)) then
            allocate(info(ng))
            info = local_status
         end if
         return
      end if
      if (present(pair_counts)) then
         if (any(shape(pair_counts) /= shape(correlations))) then
            local_status = 2
            if (present(info)) then
               allocate(info(ng))
               info = local_status
            end if
            return
         end if
      end if

      do grp = 1, ng
         if (present(pair_counts)) then
            mean_count = sum(real(pair_counts(:, :, grp), dp)) / real(size(pair_counts, 1) * size(pair_counts, 2), dp)
            if (mean_count < real(minimum_n, dp)) then
               local_status(grp) = 3
               cycle
            end if
         end if
         call score_overlap_summary(keys, correlations(:, :, grp), adj, a, av, g, item_cor, local_info)
         if (local_info /= 0) local_status(grp) = 10 + local_info
         adjusted_correlations(:, :, grp) = adj
         alpha(:, grp) = a
         average_r(:, grp) = av
         g6(:, grp) = g
         raw_cov = matmul(transpose(keys), matmul(correlations(:, :, grp), keys))
         scale_variance(:, grp) = diagonal_values(raw_cov)
      end do
      if (present(info)) then
         allocate(info(ng))
         info = local_status
      end if
   end subroutine score_by_correlations

   subroutine mediate_from_covariance(covariance, x_indices, mediator_indices, y_indices, a_paths, b_paths, &
      total_paths, direct_paths, indirect_product, indirect_difference, info)
      real(dp), intent(in) :: covariance(:, :) !! Full covariance or correlation matrix for all variables.
      integer, intent(in) :: x_indices(:) !! One-based independent-variable indices.
      integer, intent(in) :: mediator_indices(:) !! One-based mediator-variable indices.
      integer, intent(in) :: y_indices(:) !! One-based dependent-variable indices.
      real(dp), allocatable, intent(out) :: a_paths(:, :) !! X-by-mediator regression coefficients.
      real(dp), allocatable, intent(out) :: b_paths(:, :) !! Mediator-by-Y coefficients controlling for X.
      real(dp), allocatable, intent(out) :: total_paths(:, :) !! X-by-Y total-effect regression coefficients.
      real(dp), allocatable, intent(out) :: direct_paths(:, :) !! X-by-Y direct effects controlling for mediators.
      real(dp), allocatable, intent(out) :: indirect_product(:, :) !! Sum of A-times-B indirect effects for each X and Y.
      real(dp), allocatable, intent(out) :: indirect_difference(:, :) !! Total minus direct effects for each X and Y.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid indices or singular-block failures.
      real(dp), allocatable :: cxx(:, :), cxm(:, :), cxy(:, :), czz(:, :), czy(:, :), inv(:, :), beta_joint(:, :)
      integer, allocatable :: z_indices(:)
      integer :: nx, nm, ny, nz, local_info

      nx = size(x_indices)
      nm = size(mediator_indices)
      ny = size(y_indices)
      info = 0
      if (size(covariance, 1) /= size(covariance, 2) .or. nx < 1 .or. nm < 1 .or. ny < 1) info = 1
      if (info == 0) then
         if (.not. indices_valid(x_indices, size(covariance, 1))) info = 2
         if (.not. indices_valid(mediator_indices, size(covariance, 1))) info = 2
         if (.not. indices_valid(y_indices, size(covariance, 1))) info = 2
      end if
      if (info /= 0) then
         allocate(a_paths(0, 0), b_paths(0, 0), total_paths(0, 0), direct_paths(0, 0))
         allocate(indirect_product(0, 0), indirect_difference(0, 0))
         return
      end if

      call take_block(covariance, x_indices, x_indices, cxx)
      call take_block(covariance, x_indices, mediator_indices, cxm)
      call take_block(covariance, x_indices, y_indices, cxy)
      call pseudoinverse(cxx, inv, info=local_info)
      if (local_info /= 0) then
         info = 10 + local_info
         allocate(a_paths(0, 0), b_paths(0, 0), total_paths(0, 0), direct_paths(0, 0))
         allocate(indirect_product(0, 0), indirect_difference(0, 0))
         return
      end if
      allocate(a_paths(nx, nm), total_paths(nx, ny))
      a_paths = matmul(inv, cxm)
      total_paths = matmul(inv, cxy)

      nz = nx + nm
      allocate(z_indices(nz))
      z_indices(1:nx) = x_indices
      z_indices(nx + 1:nz) = mediator_indices
      call take_block(covariance, z_indices, z_indices, czz)
      call take_block(covariance, z_indices, y_indices, czy)
      call pseudoinverse(czz, inv, info=local_info)
      if (local_info /= 0) then
         info = 20 + local_info
         allocate(b_paths(0, 0), direct_paths(0, 0), indirect_product(0, 0), indirect_difference(0, 0))
         return
      end if
      allocate(beta_joint(nz, ny))
      beta_joint = matmul(inv, czy)
      allocate(direct_paths(nx, ny), b_paths(nm, ny), indirect_product(nx, ny), indirect_difference(nx, ny))
      direct_paths = beta_joint(1:nx, :)
      b_paths = beta_joint(nx + 1:nz, :)
      indirect_product = matmul(a_paths, b_paths)
      indirect_difference = total_paths - direct_paths
   end subroutine mediate_from_covariance

   subroutine factor_analysis_sapa(r, n_factors, n_obs, n_iter, reference_loadings, mean_loadings, sd_loadings, info)
      real(dp), intent(in) :: r(:, :) !! Population or observed correlation matrix used for SAPA replication.
      integer, intent(in) :: n_factors !! Number of MINRES factors fitted in every replication.
      integer, intent(in) :: n_obs !! Number of simulated observations per replication.
      integer, intent(in) :: n_iter !! Number of Monte Carlo replications used for loading uncertainty.
      real(dp), allocatable, intent(out) :: reference_loadings(:, :) !! Promax-rotated loadings fitted to the supplied matrix.
      real(dp), allocatable, intent(out) :: mean_loadings(:, :) !! Mean target-aligned loading over replications.
      real(dp), allocatable, intent(out) :: sd_loadings(:, :) !! Sample standard deviation of aligned replicated loadings.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions or numerical failures.
      real(dp), allocatable :: uniqueness(:), base(:, :), transform(:, :), phi(:, :), sample(:, :), rs(: ,:)
      real(dp), allocatable :: load(:, :), rotated(:, :), aligned(:, :), rotation(:, :), sums(:, :), sumsq(:, :)
      real(dp) :: objective
      integer :: iter, local_info, p, k

      p = size(r, 1)
      k = min(max(1, n_factors), p)
      info = 0
      if (size(r, 2) /= p .or. p < 2 .or. n_obs < 2 .or. n_iter < 1) then
         info = 1
         allocate(reference_loadings(0, 0), mean_loadings(0, 0), sd_loadings(0, 0))
         return
      end if
      call minimum_residual_factor(r, k, base, uniqueness, objective, local_info)
      if (local_info /= 0) then
         info = 10 + local_info
         allocate(reference_loadings(0, 0), mean_loadings(0, 0), sd_loadings(0, 0))
         return
      end if
      allocate(reference_loadings(p, k))
      if (k > 1) then
         allocate(transform(k, k), phi(k, k))
         call promax_rotate(base, 4, reference_loadings, transform, phi, local_info)
         if (local_info /= 0) then
            info = 20 + local_info
            return
         end if
      else
         reference_loadings = base
      end if
      allocate(sums(p, k), sumsq(p, k))
      sums = 0.0_dp
      sumsq = 0.0_dp
      do iter = 1, n_iter
         call simulate_correlated_normal(n_obs, r, sample, local_info)
         if (local_info /= 0) then
            info = 30 + local_info
            return
         end if
         allocate(rs(p, p))
         call correlation_matrix(sample, rs)
         call minimum_residual_factor(rs, k, load, uniqueness, objective, local_info)
         if (local_info /= 0) then
            info = 40 + local_info
            return
         end if
         allocate(rotated(p, k))
         if (k > 1) then
            if (.not. allocated(transform)) allocate(transform(k, k))
            if (.not. allocated(phi)) allocate(phi(k, k))
            call promax_rotate(load, 4, rotated, transform, phi, local_info)
            if (local_info /= 0) then
               info = 50 + local_info
               return
            end if
         else
            rotated = load
         end if
         allocate(aligned(p, k), rotation(k, k))
         call procrustes_rotate(rotated, reference_loadings, aligned, rotation, local_info)
         if (local_info /= 0) then
            info = 60 + local_info
            return
         end if
         sums = sums + aligned
         sumsq = sumsq + aligned * aligned
         deallocate(rs, rotated, aligned, rotation)
      end do
      allocate(mean_loadings(p, k), sd_loadings(p, k))
      mean_loadings = sums / real(n_iter, dp)
      if (n_iter > 1) then
         sd_loadings = sqrt(max(0.0_dp, (sumsq - real(n_iter, dp) * mean_loadings**2) / real(n_iter - 1, dp)))
      else
         sd_loadings = 0.0_dp
      end if
   end subroutine factor_analysis_sapa

   subroutine factor_analysis_pooled(correlations, n_factors, reference_loadings, mean_loadings, sd_loadings, info)
      real(dp), intent(in) :: correlations(:, :, :) !! Variable correlation matrices for independent datasets/replications.
      integer, intent(in) :: n_factors !! Number of MINRES factors fitted to every correlation matrix.
      real(dp), allocatable, intent(out) :: reference_loadings(:, :) !! Promax solution from the first dataset.
      real(dp), allocatable, intent(out) :: mean_loadings(:, :) !! Mean target-aligned loading across datasets.
      real(dp), allocatable, intent(out) :: sd_loadings(:, :) !! Sample standard deviation of target-aligned loadings.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions or numerical failures.
      real(dp), allocatable :: base(:, :), uniqueness(:), transform(:, :), phi(:, :), load(:, :), rotated(:, :)
      real(dp), allocatable :: aligned(:, :), rotation(:, :), sums(:, :), sumsq(:, :)
      real(dp) :: objective
      integer :: d, k, nd, p, local_info

      p = size(correlations, 1)
      nd = size(correlations, 3)
      k = min(max(1, n_factors), p)
      info = 0
      if (size(correlations, 2) /= p .or. p < 2 .or. nd < 1) then
         info = 1
         allocate(reference_loadings(0, 0), mean_loadings(0, 0), sd_loadings(0, 0))
         return
      end if
      call minimum_residual_factor(correlations(:, :, 1), k, base, uniqueness, objective, local_info)
      if (local_info /= 0) then
         info = 10 + local_info
         allocate(reference_loadings(0, 0), mean_loadings(0, 0), sd_loadings(0, 0))
         return
      end if
      allocate(reference_loadings(p, k))
      if (k > 1) then
         allocate(transform(k, k), phi(k, k))
         call promax_rotate(base, 4, reference_loadings, transform, phi, local_info)
         if (local_info /= 0) then
            info = 20 + local_info
            return
         end if
      else
         reference_loadings = base
      end if
      allocate(sums(p, k), sumsq(p, k))
      sums = 0.0_dp
      sumsq = 0.0_dp
      do d = 1, nd
         call minimum_residual_factor(correlations(:, :, d), k, load, uniqueness, objective, local_info)
         if (local_info /= 0) then
            info = 30 + local_info
            return
         end if
         allocate(rotated(p, k))
         if (k > 1) then
            if (.not. allocated(transform)) allocate(transform(k, k))
            if (.not. allocated(phi)) allocate(phi(k, k))
            call promax_rotate(load, 4, rotated, transform, phi, local_info)
            if (local_info /= 0) then
               info = 40 + local_info
               return
            end if
         else
            rotated = load
         end if
         allocate(aligned(p, k), rotation(k, k))
         call procrustes_rotate(rotated, reference_loadings, aligned, rotation, local_info)
         if (local_info /= 0) then
            info = 50 + local_info
            return
         end if
         sums = sums + aligned
         sumsq = sumsq + aligned * aligned
         deallocate(rotated, aligned, rotation)
      end do
      allocate(mean_loadings(p, k), sd_loadings(p, k))
      mean_loadings = sums / real(nd, dp)
      if (nd > 1) then
         sd_loadings = sqrt(max(0.0_dp, (sumsq - real(nd, dp) * mean_loadings**2) / real(nd - 1, dp)))
      else
         sd_loadings = 0.0_dp
      end if
   end subroutine factor_analysis_pooled

   subroutine parallel_analysis_poly(data, n_iter, observed_pc, reference_pc, observed_fa, reference_fa, &
      n_components, n_factors, info, ncat, correction)
      real(dp), intent(in) :: data(:, :) !! Observation-by-variable ordinal/mixed data matrix.
      integer, intent(in) :: n_iter !! Number of independent marginal-resampling trials.
      real(dp), allocatable, intent(out) :: observed_pc(:) !! Descending eigenvalues of the observed mixed correlation matrix.
      real(dp), allocatable, intent(out) :: reference_pc(:) !! Mean resampled mixed-correlation eigenvalues.
      real(dp), allocatable, intent(out) :: observed_fa(:) !! Eigenvalues of the observed SMC-reduced correlation matrix.
      real(dp), allocatable, intent(out) :: reference_fa(:) !! Mean resampled SMC-reduced correlation eigenvalues.
      integer, intent(out) :: n_components !! Count of observed PC eigenvalues above the resampled mean reference.
      integer, intent(out) :: n_factors !! Count of reduced-matrix eigenvalues above the resampled mean reference.
      integer, intent(out) :: info !! Zero on success; nonzero for correlation/eigen failures.
      integer, intent(in), optional :: ncat !! Maximum category count auto-classified as polytomous; default eight.
      real(dp), intent(in), optional :: correction !! Empty-cell correction for latent ordinal correlations; default 0.5.
      real(dp), allocatable :: rho(:, :), sample(:, :), sample_rho(:, :), vectors(:, :), reduced(:, :), smc(:)
      real(dp), allocatable :: values(:), sum_pc(:), sum_fa(:)
      real(dp) :: add, u
      integer :: i, j, iter, p, n, pick, local_info, categories

      n = size(data, 1)
      p = size(data, 2)
      categories = 8
      if (present(ncat)) categories = max(2, ncat)
      add = 0.5_dp
      if (present(correction)) add = correction
      info = 0
      n_components = 0
      n_factors = 0
      if (n < 2 .or. p < 2 .or. n_iter < 1) then
         info = 1
         allocate(observed_pc(0), reference_pc(0), observed_fa(0), reference_fa(0))
         return
      end if
      call mixed_correlation_matrix(data, rho, ncat=categories, correction=add, info=local_info)
      if (local_info /= 0) then
         info = 10 + local_info
         allocate(observed_pc(0), reference_pc(0), observed_fa(0), reference_fa(0))
         return
      end if
      call symmetric_eigen(rho, observed_pc, vectors, local_info, descending=.true.)
      if (local_info /= 0) then
         info = 20 + local_info
         allocate(reference_pc(0), observed_fa(0), reference_fa(0))
         return
      end if
      allocate(smc(p))
      call squared_multiple_correlations(rho, smc, local_info)
      if (local_info /= 0) smc = 1.0_dp
      allocate(reduced(p, p))
      reduced = rho
      do j = 1, p
         reduced(j, j) = smc(j)
      end do
      call symmetric_eigen(reduced, observed_fa, vectors, local_info, descending=.true.)
      if (local_info /= 0) then
         info = 30 + local_info
         allocate(reference_pc(0), reference_fa(0))
         return
      end if
      allocate(reference_pc(p), reference_fa(p), sum_pc(p), sum_fa(p), sample(n, p))
      sum_pc = 0.0_dp
      sum_fa = 0.0_dp
      do iter = 1, n_iter
         do j = 1, p
            do i = 1, n
               call random_number(u)
               pick = min(n, int(u * real(n, dp)) + 1)
               sample(i, j) = data(pick, j)
            end do
         end do
         call mixed_correlation_matrix(sample, sample_rho, ncat=categories, correction=add, info=local_info)
         if (local_info /= 0) then
            info = 40 + local_info
            return
         end if
         call symmetric_eigen(sample_rho, values, vectors, local_info, descending=.true.)
         if (local_info /= 0) then
            info = 50 + local_info
            return
         end if
         sum_pc = sum_pc + values
         call squared_multiple_correlations(sample_rho, smc, local_info)
         if (local_info /= 0) smc = 1.0_dp
         reduced = sample_rho
         do j = 1, p
            reduced(j, j) = smc(j)
         end do
         call symmetric_eigen(reduced, values, vectors, local_info, descending=.true.)
         if (local_info /= 0) then
            info = 60 + local_info
            return
         end if
         sum_fa = sum_fa + values
      end do
      reference_pc = sum_pc / real(n_iter, dp)
      reference_fa = sum_fa / real(n_iter, dp)
      n_components = eigenvalue_factor_count(observed_pc, reference_pc)
      n_factors = eigenvalue_factor_count(observed_fa, reference_fa)
   end subroutine parallel_analysis_poly

   subroutine simulate_omega_study(n_variables, n_factors, n_cases, n_trials, general_loading, results, info)
      integer, intent(in) :: n_variables !! Number of simulated items, divisible by the number of major factors.
      integer, intent(in) :: n_factors !! Number of major group factors used by the minor-factor simulator.
      integer, intent(in) :: n_cases !! Number of observations per trial; zero evaluates population correlation matrices.
      integer, intent(in) :: n_trials !! Number of independent simulation/fitting replications.
      real(dp), intent(in) :: general_loading !! Constant general-factor loading supplied to every item.
      real(dp), allocatable, intent(out) :: results(:, :) !! Trial rows: model omega, omega-h, eig1/p, one-factor omega,
                                                          !! ICLUST beta, ECV, general-loading CV, and MINRES RMS.
      integer, intent(out) :: info !! Zero on success; nonzero when a simulation, eigensolution, or fit fails.
      real(dp), allocatable :: true_loadings(:, :), model(:, :), observed(:, :), sample_r(:, :), r(:, :), values(:)
      real(dp), allocatable :: vectors(:, :), omega_loadings(:, :), group_omega(:, :), one_loadings(:, :), uniqueness(:)
      real(dp), allocatable :: clusters(:, :), merge_similarity(:), merge_alpha(:), merge_beta(:)
      integer, allocatable :: merge_left(:), merge_right(:), merge_sign(:), merge_size(:)
      real(dp) :: omega_h, omega_total, omega_limit, alpha, lambda6, ecv, objective, model_omega, one_omega
      real(dp) :: beta_value, mean_general, sd_general, rms
      real(dp) :: gvec(1)
      integer :: trial, local_info, n_merges, i, j, n_off

      info = 0
      if (n_variables < 2 .or. n_factors < 1 .or. mod(n_variables, n_factors) /= 0 .or. n_cases < 0 .or. n_trials < 1) then
         info = 1
         allocate(results(0, 0))
         return
      end if
      allocate(results(n_trials, 8))
      results = 0.0_dp
      gvec(1) = general_loading
      do trial = 1, n_trials
         call simulate_minor_factors(n_variables, n_factors, true_loadings, model, n_cases=n_cases, observed=observed, &
            sample_correlation=sample_r, general_loadings=gvec, bipolar=.false., info=local_info)
         if (local_info /= 0) then
            info = 10 + local_info
            return
         end if
         if (n_cases > 0) then
            r = sample_r
         else
            r = model
         end if
         if (abs(sum(model)) > tiny(1.0_dp)) then
            model_omega = (real(n_variables, dp) * general_loading)**2 / sum(model)
         else
            model_omega = 0.0_dp
         end if
         call omega_from_correlation(r, n_factors, omega_h, omega_total, omega_limit, alpha, lambda6, ecv, &
            omega_loadings, group_omega, local_info)
         if (local_info /= 0) then
            info = 20 + local_info
            return
         end if
         call symmetric_eigen(r, values, vectors, local_info, descending=.true.)
         if (local_info /= 0) then
            info = 30 + local_info
            return
         end if
         call minimum_residual_factor(r, 1, one_loadings, uniqueness, objective, local_info)
         if (local_info /= 0) then
            info = 40 + local_info
            return
         end if
         if (abs(sum(r)) > tiny(1.0_dp)) then
            one_omega = sum(one_loadings(:, 1))**2 / sum(r)
         else
            one_omega = 0.0_dp
         end if
         call iclust_hierarchical(r, 1, clusters, merge_left, merge_right, merge_sign, merge_similarity, &
            merge_alpha, merge_beta, merge_size, n_merges, local_info)
         if (local_info /= 0) then
            info = 50 + local_info
            return
         end if
         beta_value = 0.0_dp
         if (n_merges > 0) beta_value = merge_beta(n_merges)
         mean_general = sum(omega_loadings(:, 1)) / real(n_variables, dp)
         sd_general = 0.0_dp
         if (n_variables > 1) then
            sd_general = sqrt(sum((omega_loadings(:, 1) - mean_general)**2) / real(n_variables - 1, dp))
         end if
         rms = 0.0_dp
         n_off = 0
         do j = 2, n_variables
            do i = 1, j - 1
               rms = rms + (r(i, j) - dot_product(one_loadings(i, :), one_loadings(j, :)))**2
               n_off = n_off + 1
            end do
         end do
         if (n_off > 0) rms = sqrt(rms / real(n_off, dp))
         results(trial, 1) = model_omega
         results(trial, 2) = omega_h
         results(trial, 3) = values(1) / real(n_variables, dp)
         results(trial, 4) = one_omega
         results(trial, 5) = beta_value
         results(trial, 6) = ecv
         if (abs(mean_general) > tiny(1.0_dp)) results(trial, 7) = sd_general / abs(mean_general)
         results(trial, 8) = rms
      end do
   end subroutine simulate_omega_study


   subroutine bass_ackward_factor(r, max_factors, loadings, phi, transitions, parent_factor, info, adjust)
      real(dp), intent(in) :: r(:, :) !! Item correlation matrix used for successive one-through-K factor solutions.
      integer, intent(in) :: max_factors !! Largest factor count in the hierarchy; must be between one and matrix order.
      real(dp), allocatable, intent(out) :: loadings(:, :, :) !! Item-by-factor-by-level padded rotated loading solutions.
      real(dp), allocatable, intent(out) :: phi(:, :, :) !! Padded factor-correlation matrices for each hierarchy level.
      real(dp), allocatable, intent(out) :: transitions(:, :, :) !! Correlations between level-K and level-(K-1) scores.
      integer, allocatable, intent(out) :: parent_factor(:, :) !! For each lower-level factor, strongest absolute parent index.
      integer, allocatable, intent(out) :: info(:) !! Per-level zero success or nonzero extraction/rotation status.
      logical, intent(in), optional :: adjust !! Standardize score covariances to correlations; default true.
      real(dp), allocatable :: rinv(:, :), raw(:, :), uniqueness(:), rotated(:, :), transform(:, :), local_phi(:, :)
      real(dp), allocatable :: weights(:, :, :), cross(:, :), var_current(:), var_previous(:)
      real(dp) :: objective, denom
      integer :: k, i, j, p, local_info, best_index
      logical :: do_adjust

      p = size(r, 1)
      allocate(info(max(0, max_factors)))
      info = 0
      if (p < 2 .or. size(r, 2) /= p .or. max_factors < 1 .or. max_factors > p) then
         if (size(info) > 0) info = 1
         allocate(loadings(0, 0, 0), phi(0, 0, 0), transitions(0, 0, 0), parent_factor(0, 0))
         return
      end if
      do_adjust = .true.
      if (present(adjust)) do_adjust = adjust
      call pseudoinverse(r, rinv, info=local_info)
      if (local_info /= 0) then
         info = 2
         allocate(loadings(0, 0, 0), phi(0, 0, 0), transitions(0, 0, 0), parent_factor(0, 0))
         return
      end if
      allocate(loadings(p, max_factors, max_factors), phi(max_factors, max_factors, max_factors))
      allocate(transitions(max_factors, max_factors, max_factors), parent_factor(max_factors, max_factors))
      allocate(weights(p, max_factors, max_factors))
      loadings = 0.0_dp
      phi = 0.0_dp
      transitions = 0.0_dp
      parent_factor = 0
      weights = 0.0_dp

      do k = 1, max_factors
         call minimum_residual_factor(r, k, raw, uniqueness, objective, local_info)
         if (local_info /= 0) then
            info(k) = 10 + local_info
            cycle
         end if
         allocate(rotated(p, k), transform(k, k), local_phi(k, k))
         if (k > 1) then
            call promax_rotate(raw, 4, rotated, transform, local_phi, local_info)
            if (local_info /= 0) then
               info(k) = 20 + local_info
               deallocate(rotated, transform, local_phi)
               cycle
            end if
         else
            rotated = raw
            transform = 1.0_dp
            local_phi = 1.0_dp
         end if
         loadings(:, 1:k, k) = rotated
         phi(1:k, 1:k, k) = local_phi
         weights(:, 1:k, k) = matmul(rinv, rotated)
         if (k == 1) transitions(1, 1, 1) = 1.0_dp
         deallocate(rotated, transform, local_phi)
      end do

      do k = 2, max_factors
         if (info(k) /= 0 .or. info(k - 1) /= 0) cycle
         allocate(cross(k, k - 1), var_current(k), var_previous(k - 1))
         cross = matmul(transpose(weights(:, 1:k, k)), matmul(r, weights(:, 1:k - 1, k - 1)))
         do i = 1, k
            var_current(i) = dot_product(weights(:, i, k), matmul(r, weights(:, i, k)))
         end do
         do j = 1, k - 1
            var_previous(j) = dot_product(weights(:, j, k - 1), matmul(r, weights(:, j, k - 1)))
         end do
         do j = 1, k - 1
            do i = 1, k
               if (do_adjust) then
                  denom = sqrt(max(tiny(1.0_dp), var_current(i) * var_previous(j)))
                  transitions(i, j, k) = cross(i, j) / denom
               else
                  transitions(i, j, k) = cross(i, j)
               end if
            end do
            best_index = 1
            do i = 2, k
               if (abs(transitions(i, j, k)) > abs(transitions(best_index, j, k))) best_index = i
            end do
            parent_factor(j, k) = best_index
         end do
         deallocate(cross, var_current, var_previous)
      end do
   end subroutine bass_ackward_factor

   subroutine factor_analysis_by_groups(pooled_r, group_r, n_factors, pooled_loadings, group_loadings, mean_loadings, &
      sd_loadings, mean_phi, info)
      real(dp), intent(in) :: pooled_r(:, :) !! Pooled item correlation matrix used as the alignment target.
      real(dp), intent(in) :: group_r(:, :, :) !! Item correlation matrix for each group.
      integer, intent(in) :: n_factors !! Number of common factors fitted to pooled and group matrices.
      real(dp), allocatable, intent(out) :: pooled_loadings(:, :) !! Promax-rotated pooled factor loading target.
      real(dp), allocatable, intent(out) :: group_loadings(:, :, :) !! Target-aligned item-by-factor loadings by group.
      real(dp), allocatable, intent(out) :: mean_loadings(:, :) !! Mean aligned loading across successful groups.
      real(dp), allocatable, intent(out) :: sd_loadings(:, :) !! Sample standard deviation across successful groups.
      real(dp), allocatable, intent(out) :: mean_phi(:, :) !! Mean aligned factor-correlation matrix across groups.
      integer, allocatable, intent(out) :: info(:) !! Per-group zero success or nonzero extraction/rotation status.
      real(dp), allocatable :: base(:, :), uniqueness(:), transform(:, :), pooled_phi(:, :), raw(:, :), rotated(:, :)
      real(dp), allocatable :: local_phi(:, :), aligned(:, :), rotation(:, :), sums(:, :), sumsq(:, :), phi_sum(:, :)
      real(dp) :: objective
      integer :: g, k, ng, p, local_info, n_good

      p = size(pooled_r, 1)
      ng = size(group_r, 3)
      k = min(max(1, n_factors), p)
      allocate(info(ng))
      info = 0
      if (p < 2 .or. size(pooled_r, 2) /= p .or. size(group_r, 1) /= p .or. size(group_r, 2) /= p .or. ng < 1) then
         info = 1
         allocate(pooled_loadings(0, 0), group_loadings(0, 0, 0), mean_loadings(0, 0), sd_loadings(0, 0), mean_phi(0, 0))
         return
      end if
      call minimum_residual_factor(pooled_r, k, base, uniqueness, objective, local_info)
      if (local_info /= 0) then
         info = 2
         allocate(pooled_loadings(0, 0), group_loadings(0, 0, 0), mean_loadings(0, 0), sd_loadings(0, 0), mean_phi(0, 0))
         return
      end if
      allocate(pooled_loadings(p, k), pooled_phi(k, k), transform(k, k))
      if (k > 1) then
         call promax_rotate(base, 4, pooled_loadings, transform, pooled_phi, local_info)
         if (local_info /= 0) then
            info = 3
            allocate(group_loadings(0, 0, 0), mean_loadings(0, 0), sd_loadings(0, 0), mean_phi(0, 0))
            return
         end if
      else
         pooled_loadings = base
         pooled_phi = 1.0_dp
      end if
      allocate(group_loadings(p, k, ng), sums(p, k), sumsq(p, k), phi_sum(k, k))
      group_loadings = 0.0_dp
      sums = 0.0_dp
      sumsq = 0.0_dp
      phi_sum = 0.0_dp
      n_good = 0
      do g = 1, ng
         call minimum_residual_factor(group_r(:, :, g), k, raw, uniqueness, objective, local_info)
         if (local_info /= 0) then
            info(g) = 10 + local_info
            cycle
         end if
         allocate(rotated(p, k), local_phi(k, k))
         if (k > 1) then
            call promax_rotate(raw, 4, rotated, transform, local_phi, local_info)
            if (local_info /= 0) then
               info(g) = 20 + local_info
               deallocate(rotated, local_phi)
               cycle
            end if
         else
            rotated = raw
            local_phi = 1.0_dp
         end if
         allocate(aligned(p, k), rotation(k, k))
         call procrustes_rotate(rotated, pooled_loadings, aligned, rotation, local_info)
         if (local_info /= 0) then
            info(g) = 30 + local_info
            deallocate(rotated, local_phi, aligned, rotation)
            cycle
         end if
         group_loadings(:, :, g) = aligned
         sums = sums + aligned
         sumsq = sumsq + aligned * aligned
         phi_sum = phi_sum + matmul(transpose(rotation), matmul(local_phi, rotation))
         n_good = n_good + 1
         deallocate(rotated, local_phi, aligned, rotation)
      end do
      allocate(mean_loadings(p, k), sd_loadings(p, k), mean_phi(k, k))
      if (n_good > 0) then
         mean_loadings = sums / real(n_good, dp)
         mean_phi = phi_sum / real(n_good, dp)
      else
         mean_loadings = 0.0_dp
         mean_phi = 0.0_dp
      end if
      if (n_good > 1) then
         sd_loadings = sqrt(max(0.0_dp, (sumsq - real(n_good, dp) * mean_loadings**2) / real(n_good - 1, dp)))
      else
         sd_loadings = 0.0_dp
      end if
   end subroutine factor_analysis_by_groups

   subroutine simulate_multi_subject(n_subjects, n_variables, n_factors, n_trials, days, values, time, factor_loadings, &
      phi, beta, individual_mean, individual_sd, sine_amplitude, cosine_amplitude, ar1, info)
      integer, intent(in) :: n_subjects !! Number of independent subjects to simulate.
      integer, intent(in) :: n_variables !! Number of observed variables generated from latent factors.
      integer, intent(in) :: n_factors !! Number of within-subject latent factors.
      integer, intent(in) :: n_trials !! Number of repeated observations per subject.
      integer, intent(in) :: days !! Number of days spanned by the equally spaced trial times.
      real(dp), allocatable, intent(out) :: values(:, :, :) !! Trial-by-(factor+variable)-by-subject simulated values.
      real(dp), allocatable, intent(out) :: time(:) !! Trial times in hours from the start of the simulation.
      real(dp), intent(in), optional :: factor_loadings(:, :) !! Variable-by-factor loadings; default simple structure 0.9.
      real(dp), intent(in), optional :: phi(:, :) !! Factor correlation/mixing matrix with unit diagonal; identity by default.
      real(dp), intent(in), optional :: beta(:, :) !! Subject-by-factor linear time slopes; zeros by default.
      real(dp), intent(in), optional :: individual_mean(:, :) !! Subject-by-factor innovation means; zeros by default.
      real(dp), intent(in), optional :: individual_sd(:, :) !! Subject-by-factor innovation standard deviations; ones by default.
      real(dp), intent(in), optional :: sine_amplitude(:, :) !! Subject-by-factor circadian sine coefficients; zeros by default.
      real(dp), intent(in), optional :: cosine_amplitude(:, :) !! Subject-by-factor circadian cosine coefficients; zeros by default.
      real(dp), intent(in), optional :: ar1(:, :) !! Subject-by-factor lagged-innovation coefficients; zeros by default.
      integer, intent(out), optional :: info !! Zero on success; nonzero for incompatible dimensions or communalities above one.
      real(dp), allocatable :: loads(:, :), local_phi(:, :), local_beta(:, :), local_mean(:, :), local_sd(:, :)
      real(dp), allocatable :: local_sine(:, :), local_cosine(:, :), local_ar1(:, :), factors(:, :), errors(:), lagged(:)
      real(dp), allocatable :: between(:), noise(:)
      real(dp) :: hours, radians, uniqueness
      integer :: s, j, k, local_info, group_size

      local_info = 0
      if (n_subjects < 1 .or. n_variables < 1 .or. n_factors < 1 .or. n_trials < 2 .or. days < 1) local_info = 1
      if (present(factor_loadings)) then
         if (size(factor_loadings, 1) /= n_variables .or. size(factor_loadings, 2) /= n_factors) local_info = 2
      else
         if (mod(n_variables, n_factors) /= 0) local_info = 2
      end if
      if (present(phi)) then
         if (size(phi, 1) /= n_factors .or. size(phi, 2) /= n_factors) local_info = 3
      end if
      if (present(beta)) then
         if (size(beta, 1) /= n_subjects .or. size(beta, 2) /= n_factors) local_info = 4
      end if
      if (present(individual_mean)) then
         if (size(individual_mean, 1) /= n_subjects .or. size(individual_mean, 2) /= n_factors) local_info = 5
      end if
      if (present(individual_sd)) then
         if (size(individual_sd, 1) /= n_subjects .or. size(individual_sd, 2) /= n_factors) local_info = 6
      end if
      if (present(sine_amplitude)) then
         if (size(sine_amplitude, 1) /= n_subjects .or. size(sine_amplitude, 2) /= n_factors) local_info = 7
      end if
      if (present(cosine_amplitude)) then
         if (size(cosine_amplitude, 1) /= n_subjects .or. size(cosine_amplitude, 2) /= n_factors) local_info = 8
      end if
      if (present(ar1)) then
         if (size(ar1, 1) /= n_subjects .or. size(ar1, 2) /= n_factors) local_info = 9
      end if
      if (local_info /= 0) then
         allocate(values(0, 0, 0), time(0))
         if (present(info)) info = local_info
         return
      end if

      allocate(loads(n_variables, n_factors), local_phi(n_factors, n_factors))
      allocate(local_beta(n_subjects, n_factors), local_mean(n_subjects, n_factors), local_sd(n_subjects, n_factors))
      allocate(local_sine(n_subjects, n_factors), local_cosine(n_subjects, n_factors), local_ar1(n_subjects, n_factors))
      loads = 0.0_dp
      if (present(factor_loadings)) then
         loads = factor_loadings
      else
         group_size = n_variables / n_factors
         do j = 1, n_factors
            loads((j - 1) * group_size + 1:j * group_size, j) = 0.9_dp
         end do
      end if
      if (maxval(sum(loads * loads, dim=2)) > 1.0_dp + 1.0e-12_dp) then
         allocate(values(0, 0, 0), time(0))
         if (present(info)) info = 10
         return
      end if
      local_phi = 0.0_dp
      do j = 1, n_factors
         local_phi(j, j) = 1.0_dp
      end do
      if (present(phi)) local_phi = phi
      local_beta = 0.0_dp
      local_mean = 0.0_dp
      local_sd = 1.0_dp
      local_sine = 0.0_dp
      local_cosine = 0.0_dp
      local_ar1 = 0.0_dp
      if (present(beta)) local_beta = beta
      if (present(individual_mean)) local_mean = individual_mean
      if (present(individual_sd)) local_sd = individual_sd
      if (present(sine_amplitude)) local_sine = sine_amplitude
      if (present(cosine_amplitude)) local_cosine = cosine_amplitude
      if (present(ar1)) local_ar1 = ar1

      allocate(values(n_trials, n_factors + n_variables, n_subjects), time(n_trials))
      allocate(factors(n_trials, n_factors), errors(n_trials), lagged(n_trials), between(n_factors), noise(n_trials))
      hours = 24.0_dp * real(days, dp) / real(n_trials, dp)
      do j = 1, n_trials
         time(j) = hours * real(j, dp)
      end do
      values = 0.0_dp
      do s = 1, n_subjects
         call fill_standard_normal_local(between)
         do j = 1, n_factors
            call fill_standard_normal_local(errors)
            errors = local_mean(s, j) + local_sd(s, j) * errors
            lagged(1) = 0.0_dp
            lagged(2:) = errors(:n_trials - 1)
            do k = 1, n_trials
               radians = time(k) * acos(-1.0_dp) / 12.0_dp
               factors(k, j) = between(j) + local_beta(s, j) * time(k) / real(n_trials, dp) + &
                  sin(radians) * local_sine(s, j) + cos(radians) * local_cosine(s, j) + errors(k) + &
                  local_ar1(s, j) * lagged(k)
            end do
         end do
         factors = matmul(factors, local_phi)
         values(:, 1:n_factors, s) = factors
         do j = 1, n_variables
            uniqueness = sqrt(max(0.0_dp, 1.0_dp - sum(loads(j, :)**2)))
            call fill_standard_normal_local(noise)
            values(:, n_factors + j, s) = matmul(factors, loads(j, :)) + uniqueness * noise
         end do
      end do
      if (present(info)) info = 0
   end subroutine simulate_multi_subject

   pure subroutine best_scales_from_correlation(r, criteria, n_item, cut, keys, validity, adjusted_cut, info, overlap)
      real(dp), intent(in) :: r(:, :) !! Correlation matrix containing candidate items and criterion variables.
      integer, intent(in) :: criteria(:) !! One-based criterion-variable indices, one scale constructed per criterion.
      integer, intent(in) :: n_item !! Maximum target number of keyed items per criterion before ties.
      real(dp), intent(in) :: cut !! Minimum absolute item-criterion correlation required for inclusion.
      real(dp), allocatable, intent(out) :: keys(:, :) !! Signed unit-weight item-by-criterion key matrix.
      real(dp), allocatable, intent(out) :: validity(:) !! Correlation of each unit-weight composite with its criterion.
      real(dp), allocatable, intent(out) :: adjusted_cut(:) !! Effective per-criterion threshold after enforcing n_item.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid dimensions or empty selections.
      logical, intent(in), optional :: overlap !! Allow criterion variables as candidate items; default false.
      real(dp), allocatable :: work(:, :), abs_values(:), sorted(:), key(:)
      real(dp) :: covariance, variance, threshold
      integer :: i, j, p, ny, target, selected, local_info
      logical :: allow_overlap

      p = size(r, 1)
      ny = size(criteria)
      local_info = 0
      if (p < 2 .or. size(r, 2) /= p .or. ny < 1 .or. n_item < 1) local_info = 1
      if (local_info == 0) then
         if (.not. indices_valid(criteria, p)) local_info = 2
      end if
      if (local_info /= 0) then
         allocate(keys(0, 0), validity(0), adjusted_cut(0))
         if (present(info)) info = local_info
         return
      end if
      allow_overlap = .false.
      if (present(overlap)) allow_overlap = overlap
      allocate(work(p, p))
      work = r
      if (.not. allow_overlap) then
         do j = 1, ny
            do i = 1, ny
               work(criteria(i), criteria(j)) = 0.0_dp
            end do
         end do
      end if
      allocate(keys(p, ny), validity(ny), adjusted_cut(ny), abs_values(p), sorted(p), key(p))
      keys = 0.0_dp
      validity = 0.0_dp
      do j = 1, ny
         abs_values = abs(work(:, criteria(j)))
         sorted = abs_values
         call sort_descending_local(sorted)
         target = min(n_item, p)
         threshold = max(max(0.0_dp, cut), sorted(target))
         adjusted_cut(j) = threshold
         selected = 0
         do i = 1, p
            if (abs_values(i) >= threshold .and. abs_values(i) > 0.0_dp) then
               keys(i, j) = sign(1.0_dp, work(i, criteria(j)))
               selected = selected + 1
            end if
         end do
         if (selected == 0) then
            local_info = max(local_info, 3)
            cycle
         end if
         key = keys(:, j)
         covariance = dot_product(key, r(:, criteria(j)))
         variance = dot_product(key, matmul(r, key))
         if (variance > tiny(1.0_dp)) validity(j) = covariance / sqrt(variance)
      end do
      if (present(info)) info = local_info
   end subroutine best_scales_from_correlation


   subroutine exploratory_sem(r, x_indices, y_indices, n_factors_x, n_factors_y, loadings, structure, phi, residual, &
      rms, fit_off, info)
      real(dp), intent(in) :: r(:, :) !! Correlation matrix containing the X and Y variable sets.
      integer, intent(in) :: x_indices(:) !! One-based observed-variable indices defining the first ESEM set.
      integer, intent(in) :: y_indices(:) !! One-based observed-variable indices defining the second ESEM set.
      integer, intent(in) :: n_factors_x !! Number of exploratory factors extracted from the X-set correlation block.
      integer, intent(in) :: n_factors_y !! Number of exploratory factors extracted from the Y-set correlation block.
      real(dp), allocatable, intent(out) :: loadings(:, :) !! Extended pattern coefficients for X and Y factor sets.
      real(dp), allocatable, intent(out) :: structure(:, :) !! Extended structure coefficients for both factor sets.
      real(dp), allocatable, intent(out) :: phi(:, :) !! Correlations among all extended X and Y factors.
      real(dp), allocatable, intent(out) :: residual(:, :) !! Correlation residual matrix from the ESEM reproduced model.
      real(dp), intent(out) :: rms !! Root mean square off-diagonal residual correlation.
      real(dp), intent(out) :: fit_off !! Proportion of off-diagonal squared correlation reproduced by the model.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid indices or numerical failures.
      real(dp), allocatable :: rxx(:, :), ryy(:, :), rxall(:, :), ryall(:, :), rawx(:, :), rawy(:, :)
      real(dp), allocatable :: ux(:), uy(:), tx(:, :), ty(:, :), phix(:, :), phiy(:, :), extx(:, :), exty(:, :)
      real(dp), allocatable :: sx(:, :), sy(:, :), rinv(:, :), phiinv(:, :), model(:, :)
      real(dp) :: objx, objy, residual_ss, observed_ss
      integer :: p, kx, ky, i, j, n_off, local_info

      p = size(r, 1)
      kx = min(max(1, n_factors_x), size(x_indices))
      ky = min(max(1, n_factors_y), size(y_indices))
      info = 0
      rms = 0.0_dp
      fit_off = 0.0_dp
      if (p < 2 .or. size(r, 2) /= p .or. size(x_indices) < 2 .or. size(y_indices) < 2) info = 1
      if (info == 0) then
         if (.not. indices_valid(x_indices, p)) info = 2
         if (.not. indices_valid(y_indices, p)) info = 2
      end if
      if (info /= 0) then
         allocate(loadings(0, 0), structure(0, 0), phi(0, 0), residual(0, 0))
         return
      end if

      call take_block(r, x_indices, x_indices, rxx)
      call take_block(r, y_indices, y_indices, ryy)
      allocate(rxall(p, size(x_indices)), ryall(p, size(y_indices)))
      do j = 1, size(x_indices)
         rxall(:, j) = r(:, x_indices(j))
      end do
      do j = 1, size(y_indices)
         ryall(:, j) = r(:, y_indices(j))
      end do

      call minimum_residual_factor(rxx, kx, rawx, ux, objx, local_info)
      if (local_info /= 0) then
         info = 10 + local_info
         allocate(loadings(0, 0), structure(0, 0), phi(0, 0), residual(0, 0))
         return
      end if
      allocate(tx(kx, kx), phix(kx, kx))
      if (kx > 1) then
         allocate(extx(size(x_indices), kx))
         call promax_rotate(rawx, 4, extx, tx, phix, local_info)
         if (local_info /= 0) then
            info = 20 + local_info
            allocate(loadings(0, 0), structure(0, 0), phi(0, 0), residual(0, 0))
            return
         end if
         rawx = extx
         deallocate(extx)
      else
         phix = 1.0_dp
      end if
      call factor_extension_loadings(rxall, rxx, rawx, extx, local_info)
      if (local_info /= 0) then
         info = 30 + local_info
         allocate(loadings(0, 0), structure(0, 0), phi(0, 0), residual(0, 0))
         return
      end if
      sx = matmul(extx, phix)

      call minimum_residual_factor(ryy, ky, rawy, uy, objy, local_info)
      if (local_info /= 0) then
         info = 40 + local_info
         allocate(loadings(0, 0), structure(0, 0), phi(0, 0), residual(0, 0))
         return
      end if
      allocate(ty(ky, ky), phiy(ky, ky))
      if (ky > 1) then
         allocate(exty(size(y_indices), ky))
         call promax_rotate(rawy, 4, exty, ty, phiy, local_info)
         if (local_info /= 0) then
            info = 50 + local_info
            allocate(loadings(0, 0), structure(0, 0), phi(0, 0), residual(0, 0))
            return
         end if
         rawy = exty
         deallocate(exty)
      else
         phiy = 1.0_dp
      end if
      call factor_extension_loadings(ryall, ryy, rawy, exty, local_info)
      if (local_info /= 0) then
         info = 60 + local_info
         allocate(loadings(0, 0), structure(0, 0), phi(0, 0), residual(0, 0))
         return
      end if
      sy = matmul(exty, phiy)

      allocate(loadings(p, kx + ky), structure(p, kx + ky))
      loadings(:, 1:kx) = extx
      loadings(:, kx + 1:kx + ky) = exty
      structure(:, 1:kx) = sx
      structure(:, kx + 1:kx + ky) = sy
      call pseudoinverse(r, rinv, info=local_info)
      if (local_info /= 0) then
         info = 70 + local_info
         allocate(phi(0, 0), residual(0, 0))
         return
      end if
      allocate(phi(kx + ky, kx + ky))
      phi = matmul(transpose(structure), matmul(rinv, structure))
      phi = 0.5_dp * (phi + transpose(phi))
      do i = 1, size(phi, 1)
         phi(i, i) = 1.0_dp
      end do
      call pseudoinverse(phi, phiinv, info=local_info)
      if (local_info /= 0) then
         info = 80 + local_info
         allocate(residual(0, 0))
         return
      end if
      allocate(model(p, p), residual(p, p))
      model = matmul(matmul(structure, phiinv), transpose(structure))
      residual = r - model
      residual_ss = 0.0_dp
      observed_ss = 0.0_dp
      n_off = 0
      do j = 2, p
         do i = 1, j - 1
            residual_ss = residual_ss + residual(i, j)**2 + residual(j, i)**2
            observed_ss = observed_ss + r(i, j)**2 + r(j, i)**2
            n_off = n_off + 2
         end do
      end do
      if (n_off > 0) rms = sqrt(residual_ss / real(n_off, dp))
      if (observed_ss > tiny(1.0_dp)) fit_off = 1.0_dp - residual_ss / observed_ss
   end subroutine exploratory_sem

   subroutine correlation_fiml(data, means, covariance, correlation, objective, iterations, info, max_iter, tolerance)
      real(dp), intent(in) :: data(:, :) !! Observation-by-variable matrix; IEEE NaNs mark values missing at random.
      real(dp), allocatable, intent(out) :: means(:) !! Maximum-likelihood variable means from the normal FIML model.
      real(dp), allocatable, intent(out) :: covariance(:, :) !! Maximum-likelihood covariance matrix using divisor N.
      real(dp), allocatable, intent(out) :: correlation(:, :) !! Correlation matrix standardized from the FIML covariance.
      real(dp), intent(out) :: objective !! Average observed-data log-determinant plus Mahalanobis term, constants omitted.
      integer, intent(out) :: iterations !! Number of EM iterations performed.
      integer, intent(out) :: info !! Zero on convergence; one invalid data; two all-missing column; positive numerical status.
      integer, intent(in), optional :: max_iter !! Maximum EM iterations; default 300.
      real(dp), intent(in), optional :: tolerance !! Maximum parameter-change convergence threshold; default 1e-8.
      real(dp), allocatable :: imputed(:, :), new_mean(:), new_cov(:, :), sum_x(:), sum_xx(:, :), ex(:)
      real(dp), allocatable :: soo(:, :), smo(:, :), smm(:, :), inv(:, :), delta(:), cm(:), vmm(:, :)
      real(dp), allocatable :: evals(:), evecs(:, :), work(:, :)
      integer, allocatable :: obs(:), mis(:)
      real(dp) :: tol, change, logdet, quad, sd_i, sd_j
      integer :: n, p, i, j, row, no, nm, io, im, local_info, limit
      logical :: converged

      n = size(data, 1)
      p = size(data, 2)
      limit = 300
      if (present(max_iter)) limit = max(1, max_iter)
      tol = 1.0e-8_dp
      if (present(tolerance)) tol = max(1.0e-12_dp, tolerance)
      info = 0
      iterations = 0
      objective = huge(1.0_dp)
      if (n < 2 .or. p < 1) then
         info = 1
         allocate(means(0), covariance(0, 0), correlation(0, 0))
         return
      end if
      allocate(means(p), covariance(p, p), correlation(p, p), imputed(n, p))
      do j = 1, p
         call finite_mean_local(data(:, j), means(j), local_info)
         if (local_info /= 0) then
            info = 2
            covariance = 0.0_dp
            correlation = 0.0_dp
            return
         end if
         do i = 1, n
            if (ieee_is_nan(data(i, j))) then
               imputed(i, j) = means(j)
            else
               imputed(i, j) = data(i, j)
            end if
         end do
      end do
      covariance = 0.0_dp
      do row = 1, n
         covariance = covariance + outer_product_local(imputed(row, :) - means, imputed(row, :) - means)
      end do
      covariance = covariance / real(n, dp)
      do j = 1, p
         covariance(j, j) = max(covariance(j, j), 1.0e-8_dp)
      end do

      allocate(new_mean(p), new_cov(p, p), sum_x(p), sum_xx(p, p), ex(p), obs(p), mis(p))
      converged = .false.
      do iterations = 1, limit
         sum_x = 0.0_dp
         sum_xx = 0.0_dp
         do row = 1, n
            no = 0
            nm = 0
            do j = 1, p
               if (ieee_is_nan(data(row, j))) then
                  nm = nm + 1
                  mis(nm) = j
               else
                  no = no + 1
                  obs(no) = j
               end if
            end do
            ex = means
            if (no == 0) then
               sum_x = sum_x + ex
               sum_xx = sum_xx + outer_product_local(ex, ex) + covariance
            else if (nm == 0) then
               ex = data(row, :)
               sum_x = sum_x + ex
               sum_xx = sum_xx + outer_product_local(ex, ex)
            else
               allocate(soo(no, no), smo(nm, no), smm(nm, nm), delta(no), cm(nm), vmm(nm, nm))
               do j = 1, no
                  delta(j) = data(row, obs(j)) - means(obs(j))
                  do i = 1, no
                     soo(i, j) = covariance(obs(i), obs(j))
                  end do
               end do
               do j = 1, no
                  do i = 1, nm
                     smo(i, j) = covariance(mis(i), obs(j))
                  end do
               end do
               do j = 1, nm
                  do i = 1, nm
                     smm(i, j) = covariance(mis(i), mis(j))
                  end do
               end do
               call pseudoinverse(soo, inv, info=local_info)
               if (local_info /= 0) then
                  info = 10 + local_info
                  return
               end if
               cm = means(mis(1:nm)) + matmul(smo, matmul(inv, delta))
               do im = 1, nm
                  ex(mis(im)) = cm(im)
               end do
               do io = 1, no
                  ex(obs(io)) = data(row, obs(io))
               end do
               vmm = smm - matmul(smo, matmul(inv, transpose(smo)))
               sum_x = sum_x + ex
               sum_xx = sum_xx + outer_product_local(ex, ex)
               do j = 1, nm
                  do i = 1, nm
                     sum_xx(mis(i), mis(j)) = sum_xx(mis(i), mis(j)) + vmm(i, j)
                  end do
               end do
               deallocate(soo, smo, smm, delta, cm, vmm)
            end if
         end do
         new_mean = sum_x / real(n, dp)
         new_cov = sum_xx / real(n, dp) - outer_product_local(new_mean, new_mean)
         new_cov = 0.5_dp * (new_cov + transpose(new_cov))
         call symmetric_eigen(new_cov, evals, evecs, local_info, descending=.true.)
         if (local_info /= 0) then
            info = 20 + local_info
            return
         end if
         evals = max(evals, 1.0e-10_dp)
         allocate(work(p, p))
         work = evecs
         do j = 1, p
            work(:, j) = work(:, j) * evals(j)
         end do
         new_cov = matmul(work, transpose(evecs))
         deallocate(work)
         change = max(maxval(abs(new_mean - means)), maxval(abs(new_cov - covariance)))
         means = new_mean
         covariance = new_cov
         if (change <= tol) then
            converged = .true.
            exit
         end if
      end do
      if (.not. converged) info = 30

      objective = 0.0_dp
      do row = 1, n
         no = 0
         do j = 1, p
            if (.not. ieee_is_nan(data(row, j))) then
               no = no + 1
               obs(no) = j
            end if
         end do
         if (no == 0) cycle
         allocate(soo(no, no), delta(no))
         do j = 1, no
            delta(j) = data(row, obs(j)) - means(obs(j))
            do i = 1, no
               soo(i, j) = covariance(obs(i), obs(j))
            end do
         end do
         call pseudoinverse(soo, inv, info=local_info)
         if (local_info /= 0) then
            info = max(info, 40 + local_info)
            return
         end if
         call symmetric_eigen(soo, evals, evecs, local_info, descending=.true.)
         if (local_info /= 0 .or. minval(evals) <= 0.0_dp) then
            info = max(info, 50 + max(1, local_info))
            return
         end if
         logdet = sum(log(evals))
         quad = dot_product(delta, matmul(inv, delta))
         objective = objective + logdet + quad
         deallocate(soo, delta)
      end do
      objective = objective / real(n, dp)

      do j = 1, p
         sd_j = sqrt(max(tiny(1.0_dp), covariance(j, j)))
         do i = 1, p
            sd_i = sqrt(max(tiny(1.0_dp), covariance(i, i)))
            correlation(i, j) = covariance(i, j) / (sd_i * sd_j)
         end do
      end do
      do i = 1, p
         correlation(i, i) = 1.0_dp
      end do
   end subroutine correlation_fiml



   subroutine glb_algebraic_admm(covariance, glb, solution, iterations, primal_residual, dual_residual, info, &
      lower_bounds, upper_bounds, rho, tolerance, max_iter)
      real(dp), intent(in) :: covariance(:, :) !! Symmetric item covariance matrix used in the algebraic GLB problem.
      real(dp), intent(out) :: glb !! Greatest-lower-bound reliability estimate from the converged common diagonal.
      real(dp), allocatable, intent(out) :: solution(:) !! Optimized common-score diagonal corresponding to Rcsdp `y`.
      integer, intent(out) :: iterations !! ADMM iterations performed.
      real(dp), intent(out) :: primal_residual !! Frobenius norm of the final equality residual between split matrices.
      real(dp), intent(out) :: dual_residual !! Scaled Frobenius norm of the final PSD-projection change.
      integer, intent(out) :: info !! Zero on convergence; nonzero for invalid input, eigensolver failure, or nonconvergence.
      real(dp), intent(in), optional :: lower_bounds(:) !! Per-item lower bounds on the common diagonal; zero by default.
      real(dp), intent(in), optional :: upper_bounds(:) !! Per-item upper bounds; covariance diagonal by default.
      real(dp), intent(in), optional :: rho !! Positive ADMM penalty; scale-adaptive default is inverse mean variance.
      real(dp), intent(in), optional :: tolerance !! Relative primal/dual convergence tolerance; default 1e-8.
      integer, intent(in), optional :: max_iter !! Maximum ADMM iterations; default 10000.
      real(dp), allocatable :: lower(:), upper(:), x(:, :), z(:, :), z_old(:, :), scaled_dual(:, :), target(:, :)
      real(dp) :: penalty, tol, total_covariance, scale, norm_x, norm_z, relative_primal, relative_dual
      integer :: p, i, j, limit, local_info

      p = size(covariance, 1)
      glb = 0.0_dp
      iterations = 0
      primal_residual = huge(1.0_dp)
      dual_residual = huge(1.0_dp)
      info = 0
      if (p < 1 .or. size(covariance, 2) /= p) info = 1
      if (info == 0) then
         if (maxval(abs(covariance - transpose(covariance))) > 1.0e-10_dp * max(1.0_dp, maxval(abs(covariance)))) info = 2
      end if
      if (info /= 0) then
         allocate(solution(0))
         return
      end if
      allocate(lower(p), upper(p), solution(p))
      lower = 0.0_dp
      do i = 1, p
         upper(i) = covariance(i, i)
      end do
      if (present(lower_bounds)) then
         if (size(lower_bounds) /= p) then
            info = 3
            solution = 0.0_dp
            return
         end if
         lower = lower_bounds
      end if
      if (present(upper_bounds)) then
         if (size(upper_bounds) /= p) then
            info = 4
            solution = 0.0_dp
            return
         end if
         upper = upper_bounds
      end if
      if (any(lower > upper)) then
         info = 5
         solution = 0.0_dp
         return
      end if
      scale = sum(abs(upper)) / real(p, dp)
      penalty = 1.0_dp / max(scale, sqrt(tiny(1.0_dp)))
      if (present(rho)) penalty = rho
      if (penalty <= 0.0_dp) then
         info = 6
         solution = 0.0_dp
         return
      end if
      tol = 1.0e-8_dp
      if (present(tolerance)) tol = max(tolerance, 10.0_dp * epsilon(1.0_dp))
      limit = 10000
      if (present(max_iter)) limit = max(1, max_iter)

      allocate(x(p, p), z(p, p), z_old(p, p), scaled_dual(p, p), target(p, p))
      x = covariance
      do i = 1, p
         x(i, i) = upper(i)
      end do
      scaled_dual = 0.0_dp
      call project_psd_local(x, z, local_info)
      if (local_info /= 0) then
         info = 10 + local_info
         solution = upper
         return
      end if

      do iterations = 1, limit
         z_old = z
         target = z - scaled_dual
         do j = 1, p
            do i = 1, p
               if (i /= j) x(i, j) = covariance(i, j)
            end do
         end do
         do i = 1, p
            x(i, i) = min(upper(i), max(lower(i), target(i, i) - 1.0_dp / penalty))
         end do
         target = x + scaled_dual
         call project_psd_local(target, z, local_info)
         if (local_info /= 0) then
            info = 20 + local_info
            exit
         end if
         scaled_dual = scaled_dual + x - z
         primal_residual = frobenius_norm_local(x - z)
         dual_residual = penalty * frobenius_norm_local(z - z_old)
         norm_x = frobenius_norm_local(x)
         norm_z = frobenius_norm_local(z)
         relative_primal = primal_residual / max(1.0_dp, norm_x, norm_z)
         relative_dual = dual_residual / max(1.0_dp, penalty * frobenius_norm_local(scaled_dual))
         if (relative_primal <= tol .and. relative_dual <= tol) exit
         if (primal_residual > 10.0_dp * max(dual_residual, tiny(1.0_dp))) then
            penalty = 2.0_dp * penalty
            scaled_dual = 0.5_dp * scaled_dual
         else if (dual_residual > 10.0_dp * max(primal_residual, tiny(1.0_dp))) then
            penalty = 0.5_dp * penalty
            scaled_dual = 2.0_dp * scaled_dual
         end if
      end do
      if (iterations > limit) then
         iterations = limit
         if (info == 0) info = 30
      end if
      solution = [(x(i, i), i = 1, p)]
      total_covariance = sum(covariance)
      if (abs(total_covariance) <= tiny(1.0_dp)) then
         info = max(info, 7)
         glb = 0.0_dp
      else
         glb = (total_covariance - sum([(covariance(i, i), i = 1, p)]) + sum(solution)) / total_covariance
      end if
   end subroutine glb_algebraic_admm

   subroutine project_psd_local(a, projected, info)
      real(dp), intent(in) :: a(:, :) !! Symmetric matrix to project onto the positive-semidefinite cone.
      real(dp), intent(out) :: projected(:, :) !! Nearest PSD matrix in Frobenius norm by eigenvalue clipping.
      integer, intent(out) :: info !! Zero on success; eigensolver status otherwise.
      real(dp), allocatable :: values(:), vectors(:, :)
      integer :: i, j, k, p

      p = size(a, 1)
      call symmetric_eigen(0.5_dp * (a + transpose(a)), values, vectors, info, descending=.true.)
      if (info /= 0) then
         projected = 0.0_dp
         return
      end if
      values = max(values, 0.0_dp)
      projected = 0.0_dp
      do k = 1, p
         if (values(k) <= 0.0_dp) cycle
         do j = 1, p
            do i = 1, p
               projected(i, j) = projected(i, j) + values(k) * vectors(i, k) * vectors(j, k)
            end do
         end do
      end do
      projected = 0.5_dp * (projected + transpose(projected))
   end subroutine project_psd_local

   pure real(dp) function frobenius_norm_local(a) result(value)
      real(dp), intent(in) :: a(:, :) !! Matrix whose Frobenius norm is requested.
      value = sqrt(sum(a * a))
   end function frobenius_norm_local

   subroutine test_irt_benchmark(nvar, n_obs, theta_true, theta_estimated, parameter_summary, score_correlation, info, &
      low, high, use_normal)
      integer, intent(in) :: nvar !! Number of dichotomous items in the benchmark; must be at least two.
      integer, intent(in) :: n_obs !! Number of simulated subjects; must be at least two.
      real(dp), allocatable, intent(out) :: theta_true(:) !! Latent trait values used to generate benchmark responses.
      real(dp), allocatable, intent(out) :: theta_estimated(:) !! Two-parameter person scores using fitted item parameters.
      real(dp), allocatable, intent(out) :: parameter_summary(:, :) !! Item rows: true/fitted difficulty and discrimination.
      real(dp), intent(out) :: score_correlation !! Pearson correlation between generating and fitted person traits.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions, simulation, fitting, or scoring failure.
      real(dp), intent(in), optional :: low !! Lowest evenly spaced item difficulty; default -3.
      real(dp), intent(in), optional :: high !! Highest evenly spaced item difficulty; default 3.
      logical, intent(in), optional :: use_normal !! Use a normal-ogive generator when true; logistic by default.
      real(dp), allocatable :: difficulty_true(:), discrimination_true(:), items(:, :), fit(:)
      real(dp), allocatable :: rho(:, :), thresholds(:, :), loadings(:, :), discrimination(:, :), difficulty(:, :, :)
      real(dp), allocatable :: uniqueness(:)
      integer, allocatable :: int_items(:, :)
      real(dp) :: lo, hi, objective, denom, mean_true, mean_est
      integer :: i, local_info
      logical :: normal_link

      info = 0
      score_correlation = 0.0_dp
      if (nvar < 2 .or. n_obs < 2) then
         info = 1
         allocate(theta_true(0), theta_estimated(0), parameter_summary(0, 0))
         return
      end if
      lo = -3.0_dp
      hi = 3.0_dp
      if (present(low)) lo = low
      if (present(high)) hi = high
      normal_link = .false.
      if (present(use_normal)) normal_link = use_normal
      allocate(difficulty_true(nvar), discrimination_true(nvar))
      do i = 1, nvar
         difficulty_true(i) = lo + real(i - 1, dp) * (hi - lo) / real(nvar - 1, dp)
      end do
      discrimination_true = 1.0_dp
      if (normal_link) then
         call simulate_npn_items(n_obs, difficulty_true, discrimination_true, items, theta_true, info=local_info)
      else
         call simulate_npl_items(n_obs, difficulty_true, discrimination_true, items, theta_true, info=local_info)
      end if
      if (local_info /= 0) then
         info = 10 + local_info
         allocate(theta_estimated(0), parameter_summary(0, 0))
         return
      end if
      allocate(int_items(n_obs, nvar))
      int_items = nint(items)
      call irt_factor_analysis(int_items, 0, 1, 1, rho, thresholds, loadings, discrimination, difficulty, uniqueness, &
         objective, local_info)
      if (local_info /= 0 .or. size(difficulty, 2) < 1) then
         info = 20 + max(1, local_info)
         allocate(theta_estimated(0), parameter_summary(0, 0))
         return
      end if
      allocate(parameter_summary(nvar, 4))
      parameter_summary(:, 1) = difficulty_true
      parameter_summary(:, 2) = difficulty(:, 1, 1)
      parameter_summary(:, 3) = discrimination_true
      parameter_summary(:, 4) = discrimination(:, 1)
      if (sum(parameter_summary(:, 4)) < 0.0_dp) then
         parameter_summary(:, 2) = -parameter_summary(:, 2)
         parameter_summary(:, 4) = -parameter_summary(:, 4)
      end if
      call two_parameter_person_scores(parameter_summary(:, 2), parameter_summary(:, 4), items, theta_estimated, fit)
      mean_true = sum(theta_true) / real(n_obs, dp)
      mean_est = sum(theta_estimated) / real(n_obs, dp)
      denom = sqrt(sum((theta_true - mean_true)**2) * sum((theta_estimated - mean_est)**2))
      if (denom <= tiny(1.0_dp)) then
         info = 30
         score_correlation = 0.0_dp
      else
         score_correlation = sum((theta_true - mean_true) * (theta_estimated - mean_est)) / denom
      end if
   end subroutine test_irt_benchmark

   subroutine phi_demo_simulation(n, latent_correlation, cuts, tetrachoric, phis, yule, data, info)
      integer, intent(in) :: n !! Number of simulated normal observations; must be at least two.
      real(dp), intent(in) :: latent_correlation !! Correlation linking latent and observed normal scores; must be in [-1,1].
      real(dp), intent(in) :: cuts(:) !! Thresholds converting the observed normal score to binary item responses.
      real(dp), allocatable, intent(out) :: tetrachoric(:, :) !! Tetrachoric correlations among thresholded items.
      real(dp), allocatable, intent(out) :: phis(:, :) !! Pearson matrix with phi-to-tetrachoric values in the binary
                                                    !! upper triangle.
      real(dp), allocatable, intent(out) :: yule(:, :) !! Yule-Q correlation matrix among thresholded items.
      real(dp), allocatable, intent(out) :: data(:, :) !! Simulated latent, observed, and binary response columns.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions or correlation estimation failure.
      real(dp), allocatable :: latent(:), error(:), observed(:), binary_real(:, :), binary_phi(:, :), transformed(:, :)
      real(dp), allocatable :: ci(:, :), z_p(:, :), margins(:)
      integer, allocatable :: binary(:, :), types(:)
      real(dp) :: scale
      integer :: i, j, local_info, p

      info = 0
      p = size(cuts)
      if (n < 2 .or. p < 2 .or. abs(latent_correlation) > 1.0_dp) then
         info = 1
         allocate(tetrachoric(0, 0), phis(0, 0), yule(0, 0), data(0, 0))
         return
      end if
      allocate(latent(n), error(n), observed(n), binary(n, p), binary_real(n, p), data(n, p + 2))
      call fill_standard_normal_local(latent)
      call fill_standard_normal_local(error)
      scale = sqrt(max(0.0_dp, 1.0_dp - latent_correlation * latent_correlation))
      observed = latent_correlation * latent + scale * error
      do j = 1, p
         do i = 1, n
            binary(i, j) = merge(1, 0, observed(i) > cuts(j))
            binary_real(i, j) = real(binary(i, j), dp)
         end do
      end do
      data(:, 1) = latent
      data(:, 2) = observed
      data(:, 3:p + 2) = binary_real
      allocate(phis(p + 2, p + 2))
      call correlation_matrix(data, phis)
      allocate(types(p))
      types = 3
      call mixed_correlation_matrix(binary_real, tetrachoric, variable_types=types, correction=0.0_dp, info=local_info)
      if (local_info /= 0) info = 10 + local_info
      allocate(binary_phi(p, p), margins(p))
      binary_phi = phis(3:p + 2, 3:p + 2)
      do j = 1, p
         margins(j) = sum(binary_real(:, j)) / real(n, dp)
      end do
      call phi_to_tetrachoric_values(binary_phi, margins, transformed)
      do j = 2, p
         do i = 1, j - 1
            phis(i + 2, j + 2) = transformed(i, j)
         end do
      end do
      call yule_correlation_matrix(binary, 1.0_dp, .false., 0.05_dp, yule, ci, z_p)
   end subroutine phi_demo_simulation

   pure subroutine finite_mean_local(x, value, info)
      real(dp), intent(in) :: x(:) !! Numeric vector whose IEEE NaNs are omitted from the mean.
      real(dp), intent(out) :: value !! Mean of finite entries, or zero when every entry is missing.
      integer, intent(out) :: info !! Zero when at least one finite value exists; one for an all-missing vector.
      real(dp) :: total
      integer :: i, count_finite

      total = 0.0_dp
      count_finite = 0
      do i = 1, size(x)
         if (ieee_is_nan(x(i))) cycle
         total = total + x(i)
         count_finite = count_finite + 1
      end do
      if (count_finite > 0) then
         value = total / real(count_finite, dp)
         info = 0
      else
         value = 0.0_dp
         info = 1
      end if
   end subroutine finite_mean_local

   pure function outer_product_local(x, y) result(a)
      real(dp), intent(in) :: x(:) !! Left vector in the outer product.
      real(dp), intent(in) :: y(:) !! Right vector in the outer product.
      real(dp) :: a(size(x), size(y))
      integer :: i, j

      do j = 1, size(y)
         do i = 1, size(x)
            a(i, j) = x(i) * y(j)
         end do
      end do
   end function outer_product_local

   subroutine fill_standard_normal_local(x)
      real(dp), intent(out) :: x(:) !! Vector filled with independent standard-normal pseudo-random values.
      real(dp) :: u1, u2, radius, angle
      integer :: i

      i = 1
      do while (i <= size(x))
         call random_number(u1)
         call random_number(u2)
         u1 = max(u1, tiny(1.0_dp))
         radius = sqrt(-2.0_dp * log(u1))
         angle = 2.0_dp * acos(-1.0_dp) * u2
         x(i) = radius * cos(angle)
         if (i + 1 <= size(x)) x(i + 1) = radius * sin(angle)
         i = i + 2
      end do
   end subroutine fill_standard_normal_local

   pure subroutine sort_descending_local(x)
      real(dp), intent(inout) :: x(:) !! Numeric vector sorted in descending order in place.
      real(dp) :: temp
      integer :: i, j

      do i = 2, size(x)
         temp = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) >= temp) exit
            x(j + 1) = x(j)
            j = j - 1
         end do
         x(j + 1) = temp
      end do
   end subroutine sort_descending_local

   pure subroutine take_block(a, rows, cols, block)
      real(dp), intent(in) :: a(:, :) !! Source matrix from which a rectangular index block is extracted.
      integer, intent(in) :: rows(:) !! One-based source row indices in desired output order.
      integer, intent(in) :: cols(:) !! One-based source column indices in desired output order.
      real(dp), allocatable, intent(out) :: block(:, :) !! Extracted matrix block.
      integer :: i, j

      allocate(block(size(rows), size(cols)))
      do j = 1, size(cols)
         do i = 1, size(rows)
            block(i, j) = a(rows(i), cols(j))
         end do
      end do
   end subroutine take_block

   pure logical function indices_valid(indices, n) result(valid)
      integer, intent(in) :: indices(:) !! One-based indices to validate against a vector or matrix order.
      integer, intent(in) :: n !! Maximum legal one-based index.

      valid = .true.
      if (size(indices) < 1) valid = .false.
      if (size(indices) > 0) then
         if (minval(indices) < 1 .or. maxval(indices) > n) valid = .false.
      end if
   end function indices_valid

   pure function diagonal_values(a) result(values)
      real(dp), intent(in) :: a(:, :) !! Square matrix whose diagonal is returned.
      real(dp) :: values(min(size(a, 1), size(a, 2)))
      integer :: i

      do i = 1, size(values)
         values(i) = a(i, i)
      end do
   end function diagonal_values

end module psych_more_extended
