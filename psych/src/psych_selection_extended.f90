module psych_selection_extended
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
   use r_kinds, only : dp
   use psych_types, only : pca_result
   use psych_stats, only : correlation_matrix, finite_mean, finite_sd
   use psych_factor, only : principal_components_from_matrix
   use psych_factor_extended, only : minimum_residual_factor
   use psych_validity_extended, only : vss_complexity_fit
   use psych_simulation_algorithms, only : simulate_minor_factors
   use psych_parallel_extended, only : parallel_analysis_normal
   use psych_irt_extended, only : simulate_item_circumplex
   use psych_circular, only : circ_tests_values
   implicit none
   private
   public :: best_items_from_correlation, parcel_keys_from_correlation, cohen_d_by_groups
   public :: simulate_parallel_recovery, vss_parallel_random, simulate_circumplex_benchmark

contains

   pure subroutine best_items_from_correlation(r, criteria, cut, selected_indices, selected_correlations, counts, &
      n_item, absolute_order, info)
      real(dp), intent(in) :: r(:, :) !! Item correlation matrix; rows and columns must refer to the same items.
      integer, intent(in) :: criteria(:) !! One-based item indices used as selection criteria.
      real(dp), intent(in) :: cut !! Minimum absolute or signed correlation required for selection.
      integer, allocatable, intent(out) :: selected_indices(:, :) !! Ranked item indices by criterion; unused cells are zero.
      real(dp), allocatable, intent(out) :: selected_correlations(:, :) !! Correlations corresponding to `selected_indices`.
      integer, allocatable, intent(out) :: counts(:) !! Number of selected items for each criterion.
      integer, intent(in), optional :: n_item !! Maximum selected items in absolute-correlation mode; default 10.
      logical, intent(in), optional :: absolute_order !! Rank by absolute correlation when true; default true.
      integer, intent(out), optional :: info !! Zero on success; nonzero for a nonsquare matrix or invalid criterion index.
      real(dp), allocatable :: values(:), ranks(:)
      integer, allocatable :: indices(:)
      real(dp) :: key_value, key_rank
      integer :: c, i, j, p, cap, n_keep, key_index, local_info
      logical :: use_absolute

      p = size(r, 1)
      local_info = 0
      if (size(r, 2) /= p) local_info = 1
      do c = 1, size(criteria)
         if (criteria(c) < 1 .or. criteria(c) > p) local_info = 2
      end do
      allocate(selected_indices(p, size(criteria)), selected_correlations(p, size(criteria)), counts(size(criteria)))
      selected_indices = 0
      selected_correlations = ieee_value(0.0_dp, ieee_quiet_nan)
      counts = 0
      if (local_info /= 0) then
         if (present(info)) info = local_info
         return
      end if

      cap = 10
      if (present(n_item)) cap = max(0, n_item)
      use_absolute = .true.
      if (present(absolute_order)) use_absolute = absolute_order
      allocate(values(p), ranks(p), indices(p))
      do c = 1, size(criteria)
         n_keep = 0
         do i = 1, p
            if (i == criteria(c)) cycle
            if (ieee_is_nan(r(i, criteria(c)))) cycle
            if (use_absolute) then
               if (abs(r(i, criteria(c))) <= cut) cycle
            else
               if (r(i, criteria(c)) <= cut) cycle
            end if
            n_keep = n_keep + 1
            values(n_keep) = r(i, criteria(c))
            ranks(n_keep) = merge(abs(values(n_keep)), values(n_keep), use_absolute)
            indices(n_keep) = i
         end do
         do i = 2, n_keep
            key_value = values(i)
            key_rank = ranks(i)
            key_index = indices(i)
            j = i - 1
            do while (j >= 1)
               if (ranks(j) >= key_rank) exit
               values(j + 1) = values(j)
               ranks(j + 1) = ranks(j)
               indices(j + 1) = indices(j)
               j = j - 1
            end do
            values(j + 1) = key_value
            ranks(j + 1) = key_rank
            indices(j + 1) = key_index
         end do
         if (use_absolute) n_keep = min(n_keep, cap)
         counts(c) = n_keep
         if (n_keep > 0) then
            selected_indices(1:n_keep, c) = indices(1:n_keep)
            selected_correlations(1:n_keep, c) = values(1:n_keep)
         end if
      end do
      if (present(info)) info = 0
   end subroutine best_items_from_correlation

   pure subroutine parcel_keys_from_correlation(r, parcel_size, keys, flip, info)
      real(dp), intent(in) :: r(:, :) !! Item correlation or similarity matrix.
      integer, intent(in) :: parcel_size !! Parcel size; only sizes two and three are supported.
      real(dp), allocatable, intent(out) :: keys(:, :) !! Item-by-parcel key matrix with optional sign reversals.
      logical, intent(in), optional :: flip !! Reverse negatively aligned items within a parcel when true; default true.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid dimensions or unsupported parcel size.
      logical, allocatable :: available(:)
      logical :: do_flip
      real(dp) :: score, best_score
      integer :: i, j, k, parcel, p, n_parcels, first, second, third, local_info

      p = size(r, 1)
      local_info = 0
      if (size(r, 2) /= p .or. p < parcel_size) local_info = 1
      if (parcel_size /= 2 .and. parcel_size /= 3) local_info = 2
      n_parcels = 0
      if (local_info == 0) n_parcels = p / parcel_size
      allocate(keys(p, n_parcels))
      keys = 0.0_dp
      if (local_info /= 0) then
         if (present(info)) info = local_info
         return
      end if
      do_flip = .true.
      if (present(flip)) do_flip = flip
      allocate(available(p))
      available = .true.

      do parcel = 1, n_parcels
         first = 0
         second = 0
         third = 0
         best_score = -huge(1.0_dp)
         if (parcel_size == 2) then
            do j = 2, p
               if (.not. available(j)) cycle
               do i = 1, j - 1
                  if (.not. available(i)) cycle
                  if (ieee_is_nan(r(i, j))) cycle
                  score = abs(r(i, j))
                  if (score > best_score) then
                     best_score = score
                     first = i
                     second = j
                  end if
               end do
            end do
            if (first == 0) exit
            keys(first, parcel) = 1.0_dp
            keys(second, parcel) = 1.0_dp
            if (do_flip .and. r(first, second) < 0.0_dp) keys(second, parcel) = -1.0_dp
            available(first) = .false.
            available(second) = .false.
         else
            do k = 3, p
               if (.not. available(k)) cycle
               do j = 2, k - 1
                  if (.not. available(j)) cycle
                  do i = 1, j - 1
                     if (.not. available(i)) cycle
                     if (ieee_is_nan(r(i, j)) .or. ieee_is_nan(r(i, k)) .or. ieee_is_nan(r(j, k))) cycle
                     score = abs(r(i, j)) + abs(r(i, k)) + abs(r(j, k))
                     if (score > best_score) then
                        best_score = score
                        first = i
                        second = j
                        third = k
                     end if
                  end do
               end do
            end do
            if (first == 0) exit
            keys(first, parcel) = 1.0_dp
            keys(second, parcel) = 1.0_dp
            keys(third, parcel) = 1.0_dp
            if (do_flip .and. r(first, second) < 0.0_dp) keys(second, parcel) = -1.0_dp
            if (do_flip .and. r(first, third) < 0.0_dp) keys(third, parcel) = -1.0_dp
            available(first) = .false.
            available(second) = .false.
            available(third) = .false.
         end if
      end do
      if (present(info)) info = 0
   end subroutine parcel_keys_from_correlation

   pure subroutine cohen_d_by_groups(x, group, strata, strata_values, cohen_d, hedges_g, r_equivalent, t_value, &
      n_first, n_second, info)
      real(dp), intent(in) :: x(:, :) !! Observation-by-variable numeric data with IEEE NaNs allowed.
      integer, intent(in) :: group(:) !! Two-level comparison-group label for each observation.
      integer, intent(in) :: strata(:) !! Secondary grouping label defining independent `cohen.d.by` strata.
      integer, allocatable, intent(out) :: strata_values(:) !! Distinct stratum labels in first-observed order.
      real(dp), allocatable, intent(out) :: cohen_d(:, :) !! Stratum-by-variable Cohen d, second group minus first group.
      real(dp), allocatable, intent(out) :: hedges_g(:, :) !! Stratum-by-variable pooled-SD standardized difference.
      real(dp), allocatable, intent(out) :: r_equivalent(:, :) !! Effect-size correlation equivalent for each d value.
      real(dp), allocatable, intent(out) :: t_value(:, :) !! Equal-variance t statistic implied by Hedges g.
      integer, allocatable, intent(out) :: n_first(:, :) !! Finite observation counts for the first comparison group.
      integer, allocatable, intent(out) :: n_second(:, :) !! Finite observation counts for the second comparison group.
      integer, intent(out), optional :: info !! Zero on success; nonzero unless row counts match and exactly two groups exist.
      integer, allocatable :: unique_groups(:), unique_strata(:)
      real(dp), allocatable :: first_values(:), second_values(:)
      real(dp) :: mean1, mean2, sd1, sd2, sd_d, sd_g, p1, p2, difference
      integer :: i, j, s, n1, n2, ng, ns, local_info
      logical :: seen

      local_info = 0
      if (size(group) /= size(x, 1) .or. size(strata) /= size(x, 1)) local_info = 1
      allocate(unique_groups(size(group)), unique_strata(size(strata)))
      ng = 0
      ns = 0
      do i = 1, size(group)
         seen = .false.
         do j = 1, ng
            if (unique_groups(j) == group(i)) seen = .true.
         end do
         if (.not. seen) then
            ng = ng + 1
            unique_groups(ng) = group(i)
         end if
         seen = .false.
         do j = 1, ns
            if (unique_strata(j) == strata(i)) seen = .true.
         end do
         if (.not. seen) then
            ns = ns + 1
            unique_strata(ns) = strata(i)
         end if
      end do
      if (ng /= 2) local_info = 2
      allocate(strata_values(ns), cohen_d(ns, size(x, 2)), hedges_g(ns, size(x, 2)))
      allocate(r_equivalent(ns, size(x, 2)), t_value(ns, size(x, 2)))
      allocate(n_first(ns, size(x, 2)), n_second(ns, size(x, 2)))
      if (ns > 0) strata_values = unique_strata(1:ns)
      cohen_d = ieee_value(0.0_dp, ieee_quiet_nan)
      hedges_g = cohen_d
      r_equivalent = cohen_d
      t_value = cohen_d
      n_first = 0
      n_second = 0
      if (local_info /= 0) then
         if (present(info)) info = local_info
         return
      end if

      allocate(first_values(size(x, 1)), second_values(size(x, 1)))
      do s = 1, ns
         do j = 1, size(x, 2)
            n1 = 0
            n2 = 0
            do i = 1, size(x, 1)
               if (strata(i) /= strata_values(s) .or. ieee_is_nan(x(i, j))) cycle
               if (group(i) == unique_groups(1)) then
                  n1 = n1 + 1
                  first_values(n1) = x(i, j)
               else if (group(i) == unique_groups(2)) then
                  n2 = n2 + 1
                  second_values(n2) = x(i, j)
               end if
            end do
            n_first(s, j) = n1
            n_second(s, j) = n2
            if (n1 < 2 .or. n2 < 2) cycle
            mean1 = finite_mean(first_values(1:n1))
            mean2 = finite_mean(second_values(1:n2))
            sd1 = finite_sd(first_values(1:n1))
            sd2 = finite_sd(second_values(1:n2))
            difference = mean2 - mean1
            sd_d = sqrt(((real(n1 - 1, dp) * sd1 * sd1) + (real(n2 - 1, dp) * sd2 * sd2)) / real(n1 + n2, dp))
            sd_g = sqrt(((real(n1 - 1, dp) * sd1 * sd1) + (real(n2 - 1, dp) * sd2 * sd2)) / real(n1 + n2 - 2, dp))
            if (sd_d > tiny(1.0_dp)) cohen_d(s, j) = difference / sd_d
            if (sd_g > tiny(1.0_dp)) hedges_g(s, j) = difference / sd_g
            p1 = real(n1, dp) / real(n1 + n2, dp)
            p2 = real(n2, dp) / real(n1 + n2, dp)
            if (.not. ieee_is_nan(cohen_d(s, j)) .and. p1 * p2 > tiny(1.0_dp)) then
               r_equivalent(s, j) = cohen_d(s, j) / sqrt(cohen_d(s, j) ** 2 + 1.0_dp / (p1 * p2))
            end if
            if (.not. ieee_is_nan(hedges_g(s, j))) then
               t_value(s, j) = hedges_g(s, j) / sqrt(1.0_dp / real(n1, dp) + 1.0_dp / real(n2, dp))
            end if
         end do
      end do
      if (present(info)) info = 0
   end subroutine cohen_d_by_groups

   subroutine simulate_parallel_recovery(n_trials, n_variables, n_factors, n_subjects, n_reference_iter, results, info)
      integer, intent(in) :: n_trials !! Number of Monte Carlo trials for each design combination.
      integer, intent(in) :: n_variables(:) !! Variable counts examined by the simulation grid.
      integer, intent(in) :: n_factors(:) !! True major-factor counts examined by the simulation grid.
      integer, intent(in) :: n_subjects(:) !! Sample sizes examined by the simulation grid.
      integer, intent(in) :: n_reference_iter !! Normal-reference replicates within each parallel analysis.
      integer, allocatable, intent(out) :: results(:, :) !! Rows contain true factors, variables, trial, FA, SMC-FA, PC, subjects.
      integer, intent(out) :: info !! Zero on success; otherwise the first simulation or eigensolution failure code.
      real(dp), allocatable :: generating_loadings(:, :), model(:, :), observed(:, :), sample_r(:, :)
      real(dp), allocatable :: obs_pc(:), ref_pc_mean(:), ref_pc_q(:), obs_fa(:), ref_fa_mean(:), ref_fa_q(:)
      integer :: fi, vi, si, trial, row, nfa, nfa_smc, npc, dummy_components, local_info, n_rows

      if (n_trials < 1 .or. n_reference_iter < 1 .or. any(n_variables < 2) .or. any(n_factors < 1) .or. &
         any(n_subjects < 2)) then
         info = 1
         allocate(results(0, 7))
         return
      end if
      n_rows = n_trials * size(n_variables) * size(n_factors) * size(n_subjects)
      allocate(results(n_rows, 7))
      row = 0
      do fi = 1, size(n_factors)
         do vi = 1, size(n_variables)
            do si = 1, size(n_subjects)
               if (mod(n_variables(vi), n_factors(fi)) /= 0) then
                  info = 2
                  return
               end if
               do trial = 1, n_trials
                  call simulate_minor_factors(n_variables(vi), n_factors(fi), generating_loadings, model, &
                     n_cases=n_subjects(si), observed=observed, sample_correlation=sample_r, info=local_info)
                  if (local_info /= 0) then
                     info = 10 + local_info
                     return
                  end if
                  call parallel_analysis_normal(sample_r, n_subjects(si), n_reference_iter, 0.95_dp, obs_pc, &
                     ref_pc_mean, ref_pc_q, obs_fa, ref_fa_mean, ref_fa_q, npc, nfa, local_info, smc=.false.)
                  if (local_info /= 0) then
                     info = 20 + local_info
                     return
                  end if
                  call parallel_analysis_normal(sample_r, n_subjects(si), n_reference_iter, 0.95_dp, obs_pc, &
                     ref_pc_mean, ref_pc_q, obs_fa, ref_fa_mean, ref_fa_q, dummy_components, nfa_smc, local_info, smc=.true.)
                  if (local_info /= 0) then
                     info = 30 + local_info
                     return
                  end if
                  row = row + 1
                  results(row, :) = [n_factors(fi), n_variables(vi), trial, nfa, nfa_smc, npc, n_subjects(si)]
               end do
            end do
         end do
      end do
      info = 0
   end subroutine simulate_parallel_recovery

   subroutine vss_parallel_random(n_cases, n_variables, max_factors, fits, info)
      integer, intent(in) :: n_cases !! Number of independent normal observations to simulate.
      integer, intent(in) :: n_variables !! Number of simulated variables.
      integer, intent(in) :: max_factors !! Maximum MINRES factor count examined.
      real(dp), allocatable, intent(out) :: fits(:, :) !! Factor-count rows by VSS complexity-one/two fit columns.
      integer, intent(out) :: info !! Zero on success; otherwise invalid dimensions or factor-extraction failure.
      real(dp), allocatable :: data(:, :), r(:, :), loadings(:, :), uniqueness(:), phi(:, :), fit(:), residual_sum(:)
      real(dp) :: objective
      integer :: k, local_info, kmax

      if (n_cases < 2 .or. n_variables < 2 .or. max_factors < 1) then
         info = 1
         allocate(fits(0, 2))
         return
      end if
      kmax = min(max_factors, n_variables)
      allocate(data(n_cases, n_variables), r(n_variables, n_variables), fits(kmax, 2))
      call fill_standard_normal_local(data)
      call correlation_matrix(data, r)
      do k = 1, kmax
         allocate(loadings(n_variables, k), uniqueness(n_variables), phi(k, k))
         call minimum_residual_factor(r, k, loadings, uniqueness, objective, local_info)
         if (local_info /= 0) then
            info = 10 + local_info
            return
         end if
         phi = 0.0_dp
         phi = identity_local(k)
         call vss_complexity_fit(r, loadings, phi, 2, fit, residual_sum)
         fits(k, 1) = fit(1)
         if (size(fit) >= 2) then
            fits(k, 2) = fit(2)
         else
            fits(k, 2) = fit(1)
         end if
         deallocate(loadings, uniqueness, phi, fit, residual_sum)
      end do
      info = 0
   end subroutine vss_parallel_random

   subroutine simulate_circumplex_benchmark(sample_sizes, variable_counts, results, info)
      integer, intent(in) :: sample_sizes(:) !! Subject counts for each circumplex simulation condition.
      integer, intent(in) :: variable_counts(:) !! Item counts for each circumplex simulation condition.
      real(dp), allocatable, intent(out) :: results(:, :) !! Rows hold nvar, n, then four tests for circle, ellipse,
                                                        !! and simple structure.
      integer, intent(out) :: info !! Zero on success; nonzero when a PCA calculation fails or a count is invalid.
      real(dp), allocatable :: circle(:, :), ellipse(:, :), simple(:, :), r(:, :)
      type(pca_result) :: pc
      real(dp) :: tests(4)
      integer :: si, vi, row, n_rows

      if (any(sample_sizes < 2) .or. any(variable_counts < 4)) then
         info = 1
         allocate(results(0, 14))
         return
      end if
      n_rows = size(sample_sizes) * size(variable_counts)
      allocate(results(n_rows, 14))
      row = 0
      do si = 1, size(sample_sizes)
         do vi = 1, size(variable_counts)
            if (mod(variable_counts(vi), 4) /= 0) then
               info = 2
               return
            end if
            call simulate_item_circumplex(sample_sizes(si), variable_counts(vi), .true., 0.6_dp, 0.6_dp, 0.0_dp, &
               0.0_dp, 0.0_dp, .false., -3.0_dp, 3.0_dp, .false., 0.0_dp, circle)
            call simulate_item_circumplex(sample_sizes(si), variable_counts(vi), .true., 0.6_dp, 0.4_dp, 0.0_dp, &
               0.0_dp, 0.0_dp, .false., -3.0_dp, 3.0_dp, .false., 0.0_dp, ellipse)
            call simulate_item_circumplex(sample_sizes(si), variable_counts(vi), .false., 0.6_dp, 0.6_dp, 0.0_dp, &
               0.0_dp, 0.0_dp, .false., -3.0_dp, 3.0_dp, .false., 0.0_dp, simple)
            row = row + 1
            results(row, 1:2) = [real(variable_counts(vi), dp), real(sample_sizes(si), dp)]
            if (allocated(r)) deallocate(r)
            allocate(r(variable_counts(vi), variable_counts(vi)))
            call correlation_matrix(circle, r)
            call principal_components_from_matrix(r, 2, pc)
            if (pc%info /= 0) then
               info = 10 + pc%info
               return
            end if
            call circ_tests_values(pc%loadings, tests(1), tests(2), tests(3), tests(4))
            results(row, 3:6) = tests
            call correlation_matrix(ellipse, r)
            call principal_components_from_matrix(r, 2, pc)
            if (pc%info /= 0) then
               info = 20 + pc%info
               return
            end if
            call circ_tests_values(pc%loadings, tests(1), tests(2), tests(3), tests(4))
            results(row, 7:10) = tests
            call correlation_matrix(simple, r)
            call principal_components_from_matrix(r, 2, pc)
            if (pc%info /= 0) then
               info = 30 + pc%info
               return
            end if
            call circ_tests_values(pc%loadings, tests(1), tests(2), tests(3), tests(4))
            results(row, 11:14) = tests
         end do
      end do
      info = 0
   end subroutine simulate_circumplex_benchmark

   subroutine fill_standard_normal_local(x)
      real(dp), intent(out) :: x(:, :) !! Matrix filled with independent standard-normal pseudorandom values.
      real(dp) :: angle, radius, u1, u2
      integer :: i, j

      do j = 1, size(x, 2)
         i = 1
         do while (i <= size(x, 1))
            call random_number(u1)
            call random_number(u2)
            u1 = max(u1, tiny(1.0_dp))
            radius = sqrt(-2.0_dp * log(u1))
            angle = 2.0_dp * acos(-1.0_dp) * u2
            x(i, j) = radius * cos(angle)
            if (i + 1 <= size(x, 1)) x(i + 1, j) = radius * sin(angle)
            i = i + 2
         end do
      end do
   end subroutine fill_standard_normal_local

   pure function identity_local(n) result(identity)
      integer, intent(in) :: n !! Order of the identity matrix.
      real(dp) :: identity(n, n)
      integer :: i

      identity = 0.0_dp
      do i = 1, n
         identity(i, i) = 1.0_dp
      end do
   end function identity_local

end module psych_selection_extended
