module psych_computational_extended
   use r_kinds, only : dp
   use psych_circular, only : cosinor_fit, circadian_correlation_matrix
   use psych_correlations, only : pseudoinverse, congruence_matrix
   use psych_data_extended, only : mixed_correlation_matrix
   use psych_factor, only : promax_rotate, varimax_rotate, factor_scores_regression
   use psych_factor_extended, only : minimum_residual_factor, factor_to_irt_parameters
   use psych_irt_extended, only : irt_thresholds
   use psych_polychoric, only : polychoric_categories
   use psych_stats, only : correlation_matrix, finite_correlation, scale_columns
   use psych_types, only : cosinor_result
   use psych_parallel_extended, only : parallel_analysis_normal
   use psych_validity_extended, only : vss_complexity_fit
   implicit none
   private
   public :: multilevel_reliability_balanced, factor_analysis_poly_data
   public :: factor_analysis_multi, irt_factor_analysis, omega_from_cfa_loadings
   public :: circadian_reliability_split, cross_validation_bootstrap
   public :: factor_solution_correlation, organize_factor_loadings
   public :: parallel_select_by_keys, vss_select_by_keys

contains

   pure subroutine multilevel_reliability_balanced(values, components, reliabilities, mean_squares, info)
      real(dp), intent(in) :: values(:, :, :) !! Person-by-time-by-item balanced response array without missing cells.
      real(dp), allocatable, intent(out) :: components(:) !! Seven person/time/item variance components.
      real(dp), allocatable, intent(out) :: reliabilities(:) !! RkF, R1R, RkR, and Rc generalizability coefficients.
      real(dp), allocatable, intent(out) :: mean_squares(:) !! ANOVA mean squares in the same seven-effect order.
      integer, intent(out), optional :: info !! Zero on success; one when any factor has fewer than two levels.
      real(dp), allocatable :: mean_p(:), mean_t(:), mean_i(:)
      real(dp), allocatable :: mean_pt(:, :), mean_pi(:, :), mean_ti(:, :)
      real(dp) :: grand, ss_p, ss_t, ss_i, ss_pt, ss_pi, ss_ti, ss_e
      real(dp) :: residual, mp, mt, mi, mpt, mpi, mti, denom
      integer :: p, t, q, ip, it, iq

      p = size(values, 1)
      t = size(values, 2)
      q = size(values, 3)
      allocate(components(7), reliabilities(4), mean_squares(7))
      components = 0.0_dp
      reliabilities = 0.0_dp
      mean_squares = 0.0_dp
      if (p < 2 .or. t < 2 .or. q < 2) then
         if (present(info)) info = 1
         return
      end if
      if (present(info)) info = 0

      allocate(mean_p(p), mean_t(t), mean_i(q), mean_pt(p, t), mean_pi(p, q), mean_ti(t, q))
      grand = sum(values) / real(p * t * q, dp)
      do ip = 1, p
         mean_p(ip) = sum(values(ip, :, :)) / real(t * q, dp)
      end do
      do it = 1, t
         mean_t(it) = sum(values(:, it, :)) / real(p * q, dp)
      end do
      do iq = 1, q
         mean_i(iq) = sum(values(:, :, iq)) / real(p * t, dp)
      end do
      do ip = 1, p
         do it = 1, t
            mean_pt(ip, it) = sum(values(ip, it, :)) / real(q, dp)
         end do
      end do
      do ip = 1, p
         do iq = 1, q
            mean_pi(ip, iq) = sum(values(ip, :, iq)) / real(t, dp)
         end do
      end do
      do it = 1, t
         do iq = 1, q
            mean_ti(it, iq) = sum(values(:, it, iq)) / real(p, dp)
         end do
      end do

      ss_p = real(t * q, dp) * sum((mean_p - grand)**2)
      ss_t = real(p * q, dp) * sum((mean_t - grand)**2)
      ss_i = real(p * t, dp) * sum((mean_i - grand)**2)
      ss_pt = 0.0_dp
      do ip = 1, p
         do it = 1, t
            ss_pt = ss_pt + (mean_pt(ip, it) - mean_p(ip) - mean_t(it) + grand)**2
         end do
      end do
      ss_pt = real(q, dp) * ss_pt
      ss_pi = 0.0_dp
      do ip = 1, p
         do iq = 1, q
            ss_pi = ss_pi + (mean_pi(ip, iq) - mean_p(ip) - mean_i(iq) + grand)**2
         end do
      end do
      ss_pi = real(t, dp) * ss_pi
      ss_ti = 0.0_dp
      do it = 1, t
         do iq = 1, q
            ss_ti = ss_ti + (mean_ti(it, iq) - mean_t(it) - mean_i(iq) + grand)**2
         end do
      end do
      ss_ti = real(p, dp) * ss_ti
      ss_e = 0.0_dp
      do ip = 1, p
         do it = 1, t
            do iq = 1, q
               residual = values(ip, it, iq) - mean_pt(ip, it) - mean_pi(ip, iq) - mean_ti(it, iq) + &
                  mean_p(ip) + mean_t(it) + mean_i(iq) - grand
               ss_e = ss_e + residual * residual
            end do
         end do
      end do

      mean_squares(1) = ss_p / real(p - 1, dp)
      mean_squares(2) = ss_t / real(t - 1, dp)
      mean_squares(3) = ss_i / real(q - 1, dp)
      mean_squares(4) = ss_pt / real((p - 1) * (t - 1), dp)
      mean_squares(5) = ss_pi / real((p - 1) * (q - 1), dp)
      mean_squares(6) = ss_ti / real((t - 1) * (q - 1), dp)
      mean_squares(7) = ss_e / real((p - 1) * (t - 1) * (q - 1), dp)

      mp = mean_squares(1)
      mt = mean_squares(2)
      mi = mean_squares(3)
      mpt = mean_squares(4)
      mpi = mean_squares(5)
      mti = mean_squares(6)
      residual = mean_squares(7)
      components(1) = (mp - mpt - mpi + residual) / real(t * q, dp)
      components(2) = (mt - mpt - mti + residual) / real(p * q, dp)
      components(3) = (mi - mpi - mti + residual) / real(p * t, dp)
      components(4) = (mpt - residual) / real(q, dp)
      components(5) = (mpi - residual) / real(t, dp)
      components(6) = (mti - residual) / real(p, dp)
      components(7) = residual

      denom = components(1) + components(5) / real(q, dp) + components(7) / real(t * q, dp)
      if (abs(denom) > tiny(1.0_dp)) then
         reliabilities(1) = (components(1) + components(5) / real(q, dp)) / denom
      end if
      denom = components(1) + components(5) / real(q, dp) + components(2) + components(4) + &
         components(7) / real(q, dp)
      if (abs(denom) > tiny(1.0_dp)) then
         reliabilities(2) = (components(1) + components(5) / real(q, dp)) / denom
      end if
      denom = components(1) + components(5) / real(q, dp) + components(2) / real(t, dp) + &
         components(4) / real(t, dp) + components(7) / real(t * q, dp)
      if (abs(denom) > tiny(1.0_dp)) then
         reliabilities(3) = (components(1) + components(5) / real(q, dp)) / denom
      end if
      denom = components(4) + components(7) / real(q, dp)
      if (abs(denom) > tiny(1.0_dp)) reliabilities(4) = components(4) / denom
   end subroutine multilevel_reliability_balanced

   subroutine factor_analysis_poly_data(data, n_factors, rho, loadings, phi, uniqueness, objective, info, ncat, correction)
      real(dp), intent(in) :: data(:, :) !! Observation-by-variable numeric data for mixed latent correlations.
      integer, intent(in) :: n_factors !! Number of common factors to extract.
      real(dp), allocatable, intent(out) :: rho(:, :) !! Mixed latent/Pearson correlation matrix used for factor extraction.
      real(dp), allocatable, intent(out) :: loadings(:, :) !! Rotated MINRES loading matrix.
      real(dp), allocatable, intent(out) :: phi(:, :) !! Factor correlation matrix; identity for a one-factor solution.
      real(dp), allocatable, intent(out) :: uniqueness(:) !! MINRES uniqueness estimates.
      real(dp), intent(out) :: objective !! MINRES off-diagonal residual sum of squares.
      integer, intent(out) :: info !! Zero on success; nonzero for correlation or factor/rotation failures.
      integer, intent(in), optional :: ncat !! Maximum number of levels treated as polytomous; defaults to eight.
      real(dp), intent(in), optional :: correction !! Empty-cell correction for latent categorical correlations; defaults to 0.5.
      real(dp), allocatable :: raw_loadings(:, :), transform(:, :)
      integer :: cor_info, rot_info, max_categories
      real(dp) :: add

      max_categories = 8
      if (present(ncat)) max_categories = ncat
      add = 0.5_dp
      if (present(correction)) add = correction
      call mixed_correlation_matrix(data, rho, ncat=max_categories, correction=add, info=cor_info)
      if (cor_info /= 0) then
         allocate(loadings(size(data, 2), max(0, n_factors)))
         allocate(phi(max(0, n_factors), max(0, n_factors)))
         allocate(uniqueness(size(data, 2)))
         loadings = 0.0_dp
         phi = 0.0_dp
         uniqueness = 0.0_dp
         objective = huge(1.0_dp)
         info = cor_info
         return
      end if
      call minimum_residual_factor(rho, n_factors, raw_loadings, uniqueness, objective, info)
      if (info /= 0) return
      allocate(loadings(size(raw_loadings, 1), size(raw_loadings, 2)))
      allocate(phi(size(raw_loadings, 2), size(raw_loadings, 2)))
      if (n_factors <= 1) then
         loadings = raw_loadings
         phi = 0.0_dp
         if (n_factors == 1) phi(1, 1) = 1.0_dp
      else
         allocate(transform(n_factors, n_factors))
         call promax_rotate(raw_loadings, rotated=loadings, transformation=transform, phi=phi, info=rot_info)
         if (rot_info /= 0) then
            info = rot_info
            loadings = raw_loadings
            phi = 0.0_dp
            phi = identity_matrix(n_factors)
         end if
      end if
   end subroutine factor_analysis_poly_data

   subroutine factor_analysis_multi(r, n_first, n_second, first_loadings, first_phi, second_loadings, info)
      real(dp), intent(in) :: r(:, :) !! Item correlation matrix for the first factor stratum.
      integer, intent(in) :: n_first !! Number of first-stratum factors; must be at least two.
      integer, intent(in) :: n_second !! Number of factors extracted from the first-stratum factor-correlation matrix.
      real(dp), allocatable, intent(out) :: first_loadings(:, :) !! Promax-rotated first-stratum MINRES loadings.
      real(dp), allocatable, intent(out) :: first_phi(:, :) !! Correlations among first-stratum factors.
      real(dp), allocatable, intent(out) :: second_loadings(:, :) !! MINRES loadings of first-stratum factors on the second stratum.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions or factor/rotation failure.
      real(dp), allocatable :: raw1(:, :), uniq1(:), transform1(:, :), raw2(:, :), uniq2(:), transform2(:, :), phi2(:, :)
      real(dp) :: objective1, objective2
      integer :: rot_info

      if (n_first < 2 .or. n_second < 1 .or. n_second > n_first) then
         allocate(first_loadings(size(r, 1), max(0, n_first)))
         allocate(first_phi(max(0, n_first), max(0, n_first)))
         allocate(second_loadings(max(0, n_first), max(0, n_second)))
         first_loadings = 0.0_dp
         first_phi = 0.0_dp
         second_loadings = 0.0_dp
         info = 1
         return
      end if
      call minimum_residual_factor(r, n_first, raw1, uniq1, objective1, info)
      if (info /= 0) return
      allocate(first_loadings(size(r, 1), n_first), first_phi(n_first, n_first), transform1(n_first, n_first))
      call promax_rotate(raw1, rotated=first_loadings, transformation=transform1, phi=first_phi, info=rot_info)
      if (rot_info /= 0) then
         info = rot_info
         return
      end if
      call minimum_residual_factor(first_phi, n_second, raw2, uniq2, objective2, info)
      if (info /= 0) return
      allocate(second_loadings(n_first, n_second))
      if (n_second == 1) then
         second_loadings = raw2
      else
         allocate(transform2(n_second, n_second), phi2(n_second, n_second))
         call promax_rotate(raw2, rotated=second_loadings, transformation=transform2, phi=phi2, info=rot_info)
         if (rot_info /= 0) info = rot_info
      end if
   end subroutine factor_analysis_multi

   subroutine irt_factor_analysis(items, minimum_category, maximum_category, n_factors, rho, thresholds, loadings, &
      discrimination, difficulty, uniqueness, objective, info, correction)
      integer, intent(in) :: items(:, :) !! Subject-by-item ordered categorical responses.
      integer, intent(in) :: minimum_category !! Lowest valid response category represented in `items`.
      integer, intent(in) :: maximum_category !! Highest valid response category represented in `items`.
      integer, intent(in) :: n_factors !! Number of latent common factors.
      real(dp), allocatable, intent(out) :: rho(:, :) !! Polychoric/tetrachoric item correlation matrix.
      real(dp), allocatable, intent(out) :: thresholds(:, :) !! Item latent-normal thresholds.
      real(dp), allocatable, intent(out) :: loadings(:, :) !! Rotated MINRES factor loadings.
      real(dp), allocatable, intent(out) :: discrimination(:, :) !! IRT discrimination parameters by item and factor.
      real(dp), allocatable, intent(out) :: difficulty(:, :, :) !! IRT threshold difficulties by item, threshold, and factor.
      real(dp), allocatable, intent(out) :: uniqueness(:) !! MINRES uniqueness estimates.
      real(dp), intent(out) :: objective !! MINRES residual objective.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid inputs or factor/rotation failures.
      real(dp), intent(in), optional :: correction !! Empty-cell correction for pairwise latent correlations; defaults to 0.5.
      real(dp), allocatable :: raw_loadings(:, :), transform(:, :), phi(:, :)
      real(dp) :: add
      integer :: i, j, rot_info

      add = 0.5_dp
      if (present(correction)) add = correction
      if (size(items, 1) < 2 .or. size(items, 2) < 2 .or. n_factors < 1 .or. n_factors > size(items, 2)) then
         allocate(rho(size(items, 2), size(items, 2)), thresholds(size(items, 2), 0))
         allocate(loadings(size(items, 2), max(0, n_factors)), discrimination(size(items, 2), max(0, n_factors)))
         allocate(difficulty(size(items, 2), 0, max(0, n_factors)), uniqueness(size(items, 2)))
         rho = 0.0_dp
         loadings = 0.0_dp
         discrimination = 0.0_dp
         uniqueness = 0.0_dp
         objective = huge(1.0_dp)
         info = 1
         return
      end if
      allocate(rho(size(items, 2), size(items, 2)))
      rho = 0.0_dp
      do i = 1, size(items, 2)
         rho(i, i) = 1.0_dp
         do j = i + 1, size(items, 2)
            call polychoric_categories(items(:, i), items(:, j), rho(i, j), add)
            rho(j, i) = rho(i, j)
         end do
      end do
      call minimum_residual_factor(rho, n_factors, raw_loadings, uniqueness, objective, info)
      if (info /= 0) return
      allocate(loadings(size(raw_loadings, 1), size(raw_loadings, 2)))
      if (n_factors == 1) then
         loadings = raw_loadings
      else
         allocate(transform(n_factors, n_factors), phi(n_factors, n_factors))
         call promax_rotate(raw_loadings, rotated=loadings, transformation=transform, phi=phi, info=rot_info)
         if (rot_info /= 0) then
            info = rot_info
            loadings = raw_loadings
         end if
      end if
      call irt_thresholds(items, minimum_category, maximum_category, thresholds)
      call factor_to_irt_parameters(loadings, thresholds, discrimination, difficulty)
   end subroutine irt_factor_analysis

   subroutine omega_from_cfa_loadings(r, cfa_loadings, omega_h, omega_total, score_r2, omega_group, info)
      real(dp), intent(in) :: r(:, :) !! Item correlation matrix used to evaluate the bifactor loading solution.
      real(dp), intent(in) :: cfa_loadings(:, :) !! Item-by-factor standardized CFA loadings; first column is the general factor.
      real(dp), intent(out) :: omega_h !! Hierarchical omega attributable to the general factor.
      real(dp), intent(out) :: omega_total !! Total omega attributable to all common factors.
      real(dp), allocatable, intent(out) :: score_r2(:) !! Determinacy-like diagonal of `W' Lambda`, with `W = R^-1 Lambda`.
      real(dp), allocatable, intent(out) :: omega_group(:, :) !! Factor rows with total, general, and group-specific omega columns.
      integer, intent(out) :: info !! Zero on success; one for dimension mismatch; otherwise pseudoinverse status.
      real(dp), allocatable :: r_inv(:, :), weights(:, :)
      real(dp) :: total_variance, commonality_sum, uniqueness_sum, group_variance, group_sum, general_sum, denom
      integer, allocatable :: assignment(:)
      integer :: i, f, group_index, n_factor

      n_factor = size(cfa_loadings, 2)
      allocate(score_r2(n_factor), omega_group(n_factor, 3))
      score_r2 = 0.0_dp
      omega_group = 0.0_dp
      omega_h = 0.0_dp
      omega_total = 0.0_dp
      if (size(r, 1) /= size(r, 2) .or. size(r, 1) /= size(cfa_loadings, 1) .or. n_factor < 1) then
         info = 1
         return
      end if
      call pseudoinverse(r, r_inv, info=info)
      if (info /= 0) return
      weights = matmul(r_inv, cfa_loadings)
      do f = 1, n_factor
         score_r2(f) = dot_product(weights(:, f), cfa_loadings(:, f))
      end do
      total_variance = sum(r)
      if (abs(total_variance) <= tiny(1.0_dp)) then
         info = 2
         return
      end if
      omega_h = sum(cfa_loadings(:, 1))**2 / total_variance
      commonality_sum = sum(cfa_loadings * cfa_loadings)
      uniqueness_sum = sum([(r(i, i), i = 1, size(r, 1))]) - commonality_sum
      omega_total = (total_variance - uniqueness_sum) / total_variance
      omega_group(1, 1) = omega_total
      omega_group(1, 2) = omega_h
      if (n_factor == 1) return

      allocate(assignment(size(cfa_loadings, 1)))
      do i = 1, size(cfa_loadings, 1)
         assignment(i) = maxloc(abs(cfa_loadings(i, 2:n_factor)), dim=1) + 1
      end do
      group_sum = 0.0_dp
      do f = 2, n_factor
         group_variance = 0.0_dp
         group_sum = 0.0_dp
         general_sum = 0.0_dp
         do i = 1, size(cfa_loadings, 1)
            if (assignment(i) == f) then
               group_sum = group_sum + cfa_loadings(i, f)
               general_sum = general_sum + cfa_loadings(i, 1)
            end if
         end do
         do i = 1, size(cfa_loadings, 1)
            if (assignment(i) /= f) cycle
            do group_index = 1, size(cfa_loadings, 1)
               if (assignment(group_index) == f) group_variance = group_variance + r(i, group_index)
            end do
         end do
         denom = group_variance
         if (abs(denom) > tiny(1.0_dp)) then
            omega_group(f, 1) = (group_sum * group_sum + general_sum * general_sum) / denom
            omega_group(f, 2) = general_sum * general_sum / denom
            omega_group(f, 3) = group_sum * group_sum / denom
         end if
      end do
   end subroutine omega_from_cfa_loadings

   pure subroutine circadian_reliability_split(angle, responses, phase_reliability, fit_reliability, splits, odd_even, &
      hours, period, minimum_observations, info)
      real(dp), intent(in) :: angle(:, :) !! Subject-by-observation measurement times in radians or clock units.
      real(dp), intent(in) :: responses(:, :, :) !! Subject-by-observation-by-variable response array aligned with `angle`.
      real(dp), allocatable, intent(out) :: phase_reliability(:) !! Adjusted circular split-half phase reliability.
      real(dp), allocatable, intent(out) :: fit_reliability(:) !! Adjusted split-half reliability of cosinor fit.
      real(dp), allocatable, intent(out) :: splits(:, :, :) !! Subject-by-variable-by-four array: phase1, phase2, fit1, fit2.
      logical, intent(in), optional :: odd_even !! Use odd/even observations instead of first/second halves; default false.
      logical, intent(in), optional :: hours !! Interpret `angle` and phase output as clock units; default true.
      real(dp), intent(in), optional :: period !! Clock period used when `hours` is true; defaults to 24.
      integer, intent(in), optional :: minimum_observations !! Minimum observations required per subject; defaults to 16.
      integer, intent(out), optional :: info !! Zero on success; one for shape mismatch; two for too few observations.
      type(cosinor_result) :: fit1, fit2
      real(dp), allocatable :: phase_pair(:, :), phase_cor(:, :), fit1_values(:), fit2_values(:)
      real(dp), allocatable :: angle1(:), angle2(:), y1(:), y2(:)
      real(dp) :: r, fit_r, use_period
      integer :: s, i, j, n, n1, n2, min_obs, n_subjects, n_variables
      logical :: use_odd_even, use_hours

      n_subjects = size(responses, 1)
      n = size(responses, 2)
      n_variables = size(responses, 3)
      allocate(phase_reliability(n_variables), fit_reliability(n_variables))
      allocate(splits(n_subjects, n_variables, 4))
      phase_reliability = 0.0_dp
      fit_reliability = 0.0_dp
      splits = 0.0_dp
      if (size(angle, 1) /= n_subjects .or. size(angle, 2) /= n) then
         if (present(info)) info = 1
         return
      end if
      min_obs = 16
      if (present(minimum_observations)) min_obs = minimum_observations
      if (n < min_obs .or. n < 4) then
         if (present(info)) info = 2
         return
      end if
      if (present(info)) info = 0
      use_odd_even = .false.
      if (present(odd_even)) use_odd_even = odd_even
      use_hours = .true.
      if (present(hours)) use_hours = hours
      use_period = 24.0_dp
      if (present(period)) use_period = period

      if (use_odd_even) then
         n1 = (n + 1) / 2
         n2 = n / 2
      else
         n1 = n / 2
         n2 = n - n1
      end if
      allocate(angle1(n1), angle2(n2), y1(n1), y2(n2))
      do s = 1, n_subjects
         do j = 1, n_variables
            if (use_odd_even) then
               n1 = 0
               n2 = 0
               do i = 1, n
                  if (mod(i, 2) == 1) then
                     n1 = n1 + 1
                     angle1(n1) = angle(s, i)
                     y1(n1) = responses(s, i, j)
                  else
                     n2 = n2 + 1
                     angle2(n2) = angle(s, i)
                     y2(n2) = responses(s, i, j)
                  end if
               end do
            else
               n1 = n / 2
               n2 = n - n1
               angle1 = angle(s, 1:n1)
               y1 = responses(s, 1:n1, j)
               angle2 = angle(s, n1 + 1:n)
               y2 = responses(s, n1 + 1:n, j)
            end if
            call cosinor_fit(angle1(1:n1), y1(1:n1), fit1, hours=use_hours, period=use_period)
            call cosinor_fit(angle2(1:n2), y2(1:n2), fit2, hours=use_hours, period=use_period)
            splits(s, j, :) = [fit1%phase, fit2%phase, fit1%fit, fit2%fit]
         end do
      end do

      allocate(phase_pair(n_subjects, 2), fit1_values(n_subjects), fit2_values(n_subjects))
      do j = 1, n_variables
         phase_pair(:, 1) = splits(:, j, 1)
         phase_pair(:, 2) = splits(:, j, 2)
         call circadian_correlation_matrix(phase_pair, phase_cor, hours=use_hours, period=use_period)
         r = phase_cor(1, 2)
         fit1_values = splits(:, j, 3)
         fit2_values = splits(:, j, 4)
         fit_r = finite_correlation(fit1_values, fit2_values)
         if (abs(1.0_dp + r) > tiny(1.0_dp)) phase_reliability(j) = 2.0_dp * r / (1.0_dp + r)
         if (abs(1.0_dp + fit_r) > tiny(1.0_dp)) fit_reliability(j) = 2.0_dp * fit_r / (1.0_dp + fit_r)
      end do
   end subroutine circadian_reliability_split

   subroutine cross_validation_bootstrap(y, x, n_iter, trial_fit, mean_fit, mean_coefficients, info)
      real(dp), intent(in) :: y(:, :) !! Observation-by-criterion matrix with complete numeric values.
      real(dp), intent(in) :: x(:, :) !! Observation-by-predictor matrix aligned with `y`.
      integer, intent(in) :: n_iter !! Number of bootstrap derivation / out-of-bag validation trials.
      real(dp), allocatable, intent(out) :: trial_fit(:, :) !! Derivation and validation correlations by bootstrap trial.
      real(dp), allocatable, intent(out) :: mean_fit(:) !! Column means of derivation and cross-validation correlations.
      real(dp), allocatable, intent(out) :: mean_coefficients(:, :) !! Mean standardized regression weights.
      integer, intent(out) :: info !! Zero on success; one for incompatible dimensions; two if no trial has out-of-bag observations.
      real(dp), allocatable :: xb(:, :), xbs(:, :), yb(:, :), combined(:, :), r(:, :), rxx(:, :), rxy(:, :)
      real(dp), allocatable :: pinv(:, :), beta(:, :), coefficients(:, :, :), pred(:), u(:), means(:), sds(:), scaled(:)
      integer, allocatable :: sampled(:), counts(:), oob(:)
      logical, allocatable :: valid(:)
      real(dp) :: draw
      integer :: i, j, k, n, npred, ncrit, noob, pinv_info, good_trials

      n = size(x, 1)
      npred = size(x, 2)
      ncrit = size(y, 2)
      allocate(trial_fit(max(0, n_iter), 2 * ncrit), mean_fit(2 * ncrit), mean_coefficients(npred, ncrit))
      trial_fit = 0.0_dp
      mean_fit = 0.0_dp
      mean_coefficients = 0.0_dp
      if (size(y, 1) /= n .or. n < 3 .or. npred < 1 .or. ncrit < 1 .or. n_iter < 1) then
         info = 1
         return
      end if
      allocate(coefficients(npred, ncrit, n_iter), sampled(n), counts(n), oob(n), valid(n_iter))
      coefficients = 0.0_dp
      valid = .false.
      good_trials = 0
      do k = 1, n_iter
         counts = 0
         do i = 1, n
            call random_number(draw)
            sampled(i) = min(n, int(draw * real(n, dp)) + 1)
            counts(sampled(i)) = counts(sampled(i)) + 1
         end do
         noob = 0
         do i = 1, n
            if (counts(i) == 0) then
               noob = noob + 1
               oob(noob) = i
            end if
         end do
         if (noob < 2) cycle
         allocate(xb(n, npred), xbs(n, npred), yb(n, ncrit), combined(n, npred + ncrit), &
            r(npred + ncrit, npred + ncrit), means(npred), sds(npred), scaled(npred))
         do i = 1, n
            xb(i, :) = x(sampled(i), :)
            yb(i, :) = y(sampled(i), :)
         end do
         call scale_columns(xb, xbs, means, sds)
         combined(:, 1:npred) = xbs
         combined(:, npred + 1:) = yb
         call correlation_matrix(combined, r)
         rxx = r(1:npred, 1:npred)
         rxy = r(1:npred, npred + 1:npred + ncrit)
         call pseudoinverse(rxx, pinv, info=pinv_info)
         if (pinv_info /= 0) then
            deallocate(xb, xbs, yb, combined, r, means, sds, scaled)
            cycle
         end if
         beta = matmul(pinv, rxy)
         coefficients(:, :, k) = beta
         do j = 1, ncrit
            allocate(pred(n), u(noob))
            pred = matmul(xbs, beta(:, j))
            trial_fit(k, j) = finite_correlation(pred, yb(:, j))
            do i = 1, noob
               where (sds > tiny(1.0_dp))
                  scaled = (x(oob(i), :) - means) / sds
               elsewhere
                  scaled = 0.0_dp
               end where
               u(i) = dot_product(scaled, beta(:, j))
            end do
            trial_fit(k, ncrit + j) = finite_correlation(u, y(oob(1:noob), j))
            deallocate(pred, u)
         end do
         good_trials = good_trials + 1
         valid(k) = .true.
         deallocate(xb, xbs, yb, combined, r, means, sds, scaled)
      end do
      if (good_trials == 0) then
         info = 2
         return
      end if
      info = 0
      do j = 1, 2 * ncrit
         mean_fit(j) = sum(trial_fit(:, j), mask=valid) / real(good_trials, dp)
      end do
      do j = 1, ncrit
         do i = 1, npred
            mean_coefficients(i, j) = sum(coefficients(i, j, :), mask=valid) / real(good_trials, dp)
         end do
      end do
   end subroutine cross_validation_bootstrap


   subroutine factor_solution_correlation(r, n_first, n_second, factor_correlation, congruence, first_loadings, &
      second_loadings, first_score_variances, second_score_variances, info, adjust_first, adjust_second)
      real(dp), intent(in) :: r(:, :) !! Observed-variable correlation matrix shared by both factor solutions.
      integer, intent(in) :: n_first !! Positive MINRES factor count for the first solution; must be less than item count.
      integer, intent(in) :: n_second !! Positive MINRES factor count for the second solution; must be less than item count.
      real(dp), allocatable, intent(out) :: factor_correlation(:, :) !! Cross-correlation matrix between factor-score estimates.
      real(dp), allocatable, intent(out) :: congruence(:, :) !! Tucker congruence coefficients between the two loading solutions.
      real(dp), allocatable, intent(out) :: first_loadings(:, :) !! Rotated loadings for the first factor solution.
      real(dp), allocatable, intent(out) :: second_loadings(:, :) !! Rotated loadings for the second factor solution.
      real(dp), allocatable, intent(out) :: first_score_variances(:) !! Variances of first-solution factor-score estimates.
      real(dp), allocatable, intent(out) :: second_score_variances(:) !! Variances of second-solution factor-score estimates.
      integer, intent(out) :: info !! Zero on success; positive values mark invalid dimensions or failed factor/score steps.
      logical, intent(in), optional :: adjust_first !! Normalize rows by first-score SDs like R `adjust[1]`; default true.
      logical, intent(in), optional :: adjust_second !! Normalize columns by second-score SDs like R `adjust[2]`; default true.
      real(dp), allocatable :: raw_first(:, :), raw_second(:, :), uniqueness(:), weights_first(:, :), weights_second(:, :)
      real(dp), allocatable :: transform(:, :), phi(:, :), score_cov(:, :), cross_cov(:, :)
      real(dp) :: objective
      integer :: i, local_info, p
      logical :: do_adjust_first, do_adjust_second

      p = size(r, 1)
      do_adjust_first = .true.
      do_adjust_second = .true.
      if (present(adjust_first)) do_adjust_first = adjust_first
      if (present(adjust_second)) do_adjust_second = adjust_second
      if (p < 2 .or. size(r, 2) /= p .or. n_first < 1 .or. n_second < 1 .or. n_first >= p .or. n_second >= p) then
         info = 1
         allocate(factor_correlation(0, 0), congruence(0, 0), first_loadings(0, 0), second_loadings(0, 0))
         allocate(first_score_variances(0), second_score_variances(0))
         return
      end if

      call minimum_residual_factor(r, n_first, raw_first, uniqueness, objective, local_info)
      if (local_info /= 0) then
         info = 10 + local_info
         return
      end if
      allocate(first_loadings(p, n_first))
      if (n_first > 1) then
         allocate(transform(n_first, n_first), phi(n_first, n_first))
         call promax_rotate(raw_first, rotated=first_loadings, transformation=transform, phi=phi, info=local_info)
         if (local_info /= 0) then
            info = 20 + local_info
            return
         end if
         deallocate(transform, phi)
      else
         first_loadings = raw_first
      end if

      call minimum_residual_factor(r, n_second, raw_second, uniqueness, objective, local_info)
      if (local_info /= 0) then
         info = 30 + local_info
         return
      end if
      allocate(second_loadings(p, n_second))
      if (n_second > 1) then
         allocate(transform(n_second, n_second), phi(n_second, n_second))
         call promax_rotate(raw_second, rotated=second_loadings, transformation=transform, phi=phi, info=local_info)
         if (local_info /= 0) then
            info = 40 + local_info
            return
         end if
      else
         second_loadings = raw_second
      end if

      allocate(weights_first(p, n_first), weights_second(p, n_second))
      call factor_scores_regression(r, first_loadings, weights_first, local_info)
      if (local_info /= 0) then
         info = 50 + local_info
         return
      end if
      call factor_scores_regression(r, second_loadings, weights_second, local_info)
      if (local_info /= 0) then
         info = 60 + local_info
         return
      end if

      cross_cov = matmul(matmul(transpose(weights_first), r), weights_second)
      score_cov = matmul(matmul(transpose(weights_first), r), weights_first)
      allocate(first_score_variances(n_first))
      do i = 1, n_first
         first_score_variances(i) = score_cov(i, i)
      end do
      score_cov = matmul(matmul(transpose(weights_second), r), weights_second)
      allocate(second_score_variances(n_second))
      do i = 1, n_second
         second_score_variances(i) = score_cov(i, i)
      end do

      allocate(factor_correlation(n_first, n_second), congruence(n_first, n_second))
      factor_correlation = cross_cov
      if (do_adjust_first) then
         do i = 1, n_first
            if (first_score_variances(i) > tiny(1.0_dp)) then
               factor_correlation(i, :) = factor_correlation(i, :) / sqrt(first_score_variances(i))
            end if
         end do
      end if
      if (do_adjust_second) then
         do i = 1, n_second
            if (second_score_variances(i) > tiny(1.0_dp)) then
               factor_correlation(:, i) = factor_correlation(:, i) / sqrt(second_score_variances(i))
            end if
         end do
      end if
      call congruence_matrix(first_loadings, second_loadings, congruence)
      info = 0
   end subroutine factor_solution_correlation

   pure subroutine organize_factor_loadings(loadings, organized, column_signs, info, factor_order, item_order, flip)
      real(dp), intent(in) :: loadings(:, :) !! Numeric loading matrix to sign-normalize and optionally reorder.
      real(dp), allocatable, intent(out) :: organized(:, :) !! Reordered matrix matching R `fa.organize` for matrix input.
      real(dp), allocatable, intent(out) :: column_signs(:) !! Applied signs in final factor order; each is plus or minus one.
      integer, intent(out) :: info !! Zero on success; one for an invalid factor order and two for an invalid item order.
      integer, intent(in), optional :: factor_order(:) !! One-based permutation/subset of factor columns; default keeps input order.
      integer, intent(in), optional :: item_order(:) !! One-based permutation/subset of item rows; default keeps input order.
      logical, intent(in), optional :: flip !! Make each selected factor's loading sum nonnegative when true; default true.
      real(dp), allocatable :: signed_loadings(:, :), signs(:)
      integer, allocatable :: cols(:), rows(:)
      integer :: i, j, nf, ni
      logical :: do_flip

      nf = size(loadings, 2)
      ni = size(loadings, 1)
      do_flip = .true.
      if (present(flip)) do_flip = flip
      allocate(cols(nf), rows(ni), signs(nf), signed_loadings(ni, nf))
      cols = [(i, i = 1, nf)]
      rows = [(i, i = 1, ni)]
      if (present(factor_order)) then
         deallocate(cols)
         allocate(cols(size(factor_order)))
         cols = factor_order
         if (any(cols < 1) .or. any(cols > nf)) then
            info = 1
            allocate(organized(0, 0), column_signs(0))
            return
         end if
      end if
      if (present(item_order)) then
         deallocate(rows)
         allocate(rows(size(item_order)))
         rows = item_order
         if (any(rows < 1) .or. any(rows > ni)) then
            info = 2
            allocate(organized(0, 0), column_signs(0))
            return
         end if
      end if

      signs = 1.0_dp
      if (do_flip) then
         do j = 1, nf
            if (sum(loadings(:, j)) < 0.0_dp) signs(j) = -1.0_dp
         end do
      end if
      do j = 1, nf
         signed_loadings(:, j) = loadings(:, j) * signs(j)
      end do
      allocate(organized(size(rows), size(cols)), column_signs(size(cols)))
      do j = 1, size(cols)
         column_signs(j) = signs(cols(j))
         do i = 1, size(rows)
            organized(i, j) = signed_loadings(rows(i), cols(j))
         end do
      end do
      info = 0
   end subroutine organize_factor_loadings

   subroutine parallel_select_by_keys(data, keys, n_iter, quantile, n_factors, n_components, scale_info, smc)
      real(dp), intent(in) :: data(:, :) !! Observation-by-item numeric data matrix used to form keyed subset correlations.
      integer, intent(in) :: keys(:, :) !! Item-by-scale key matrix; any nonzero entry includes that item in the scale.
      integer, intent(in) :: n_iter !! Number of normal-reference replicates for each keyed parallel analysis.
      real(dp), intent(in) :: quantile !! Reference eigenvalue quantile in `[0,1]`, commonly 0.95.
      integer, allocatable, intent(out) :: n_factors(:) !! Parallel-analysis factor counts, one value per key column.
      integer, allocatable, intent(out) :: n_components(:) !! Parallel-analysis component counts, one value per key column.
      integer, allocatable, intent(out) :: scale_info(:) !! Zero for each successful scale; nonzero marks subset/numerical failure.
      logical, intent(in), optional :: smc !! Use SMC factor-reference diagonals when true; default false.
      real(dp), allocatable :: subset(:, :), r(:, :)
      real(dp), allocatable :: observed_pc(:), reference_pc_mean(:), reference_pc_quantile(:)
      real(dp), allocatable :: observed_fa(:), reference_fa_mean(:), reference_fa_quantile(:)
      integer :: i, j, local_info, n_selected, s
      logical :: use_smc

      s = size(keys, 2)
      allocate(n_factors(s), n_components(s), scale_info(s))
      n_factors = -1
      n_components = -1
      scale_info = 0
      use_smc = .false.
      if (present(smc)) use_smc = smc
      if (size(keys, 1) /= size(data, 2) .or. size(data, 1) < 2 .or. n_iter < 1) then
         scale_info = 1
         return
      end if
      do j = 1, s
         n_selected = count(keys(:, j) /= 0)
         if (n_selected < 2) then
            scale_info(j) = 2
            cycle
         end if
         allocate(subset(size(data, 1), n_selected), r(n_selected, n_selected))
         n_selected = 0
         do i = 1, size(data, 2)
            if (keys(i, j) /= 0) then
               n_selected = n_selected + 1
               subset(:, n_selected) = data(:, i)
            end if
         end do
         call correlation_matrix(subset, r)
         call parallel_analysis_normal(r, size(data, 1), n_iter, quantile, observed_pc, reference_pc_mean, &
            reference_pc_quantile, observed_fa, reference_fa_mean, reference_fa_quantile, n_components(j), &
            n_factors(j), local_info, smc=use_smc)
         scale_info(j) = local_info
         deallocate(subset, r)
      end do
   end subroutine parallel_select_by_keys

   subroutine vss_select_by_keys(data, keys, max_factors, vss_one, vss_two, scale_info)
      real(dp), intent(in) :: data(:, :) !! Observation-by-item numeric data matrix used to form keyed subset correlations.
      integer, intent(in) :: keys(:, :) !! Item-by-scale key matrix; any nonzero entry includes that item in the scale.
      integer, intent(in) :: max_factors !! Maximum MINRES factor count examined for each keyed subset; R defaults to eight.
      integer, allocatable, intent(out) :: vss_one(:) !! First factor count attaining maximum complexity-one VSS fit for each scale.
      integer, allocatable, intent(out) :: vss_two(:) !! First factor count attaining maximum complexity-two VSS fit for each scale.
      integer, allocatable, intent(out) :: scale_info(:) !! Zero for each successful scale; nonzero marks subset/factor failure.
      real(dp), allocatable :: subset(:, :), r(:, :), loadings(:, :), uniqueness(:), rotated(:, :)
      real(dp), allocatable :: rotation(:, :), phi(:, :), fit(:), residual_sum(:)
      real(dp) :: best_one, best_two, objective
      integer :: i, j, local_info, n_selected, nf, nf_limit, s

      s = size(keys, 2)
      allocate(vss_one(s), vss_two(s), scale_info(s))
      vss_one = 0
      vss_two = 0
      scale_info = 0
      if (size(keys, 1) /= size(data, 2) .or. size(data, 1) < 2 .or. max_factors < 1) then
         scale_info = 1
         return
      end if
      do j = 1, s
         n_selected = count(keys(:, j) /= 0)
         if (n_selected < 2) then
            scale_info(j) = 2
            cycle
         end if
         allocate(subset(size(data, 1), n_selected), r(n_selected, n_selected))
         n_selected = 0
         do i = 1, size(data, 2)
            if (keys(i, j) /= 0) then
               n_selected = n_selected + 1
               subset(:, n_selected) = data(:, i)
            end if
         end do
         call correlation_matrix(subset, r)
         nf_limit = min(max_factors, max(1, n_selected - 1))
         best_one = -huge(1.0_dp)
         best_two = -huge(1.0_dp)
         do nf = 1, nf_limit
            call minimum_residual_factor(r, nf, loadings, uniqueness, objective, local_info)
            if (local_info /= 0) then
               scale_info(j) = 10 + local_info
               exit
            end if
            allocate(rotated(n_selected, nf), phi(nf, nf))
            phi = identity_matrix(nf)
            if (nf > 1) then
               allocate(rotation(nf, nf))
               call varimax_rotate(loadings, rotated, rotation, local_info)
               deallocate(rotation)
               if (local_info /= 0) then
                  scale_info(j) = 20 + local_info
                  deallocate(rotated, phi)
                  exit
               end if
            else
               rotated = loadings
            end if
            call vss_complexity_fit(r, rotated, phi, min(2, nf), fit, residual_sum)
            if (fit(1) > best_one) then
               best_one = fit(1)
               vss_one(j) = nf
            end if
            if (nf >= 2) then
               if (fit(2) > best_two) then
                  best_two = fit(2)
                  vss_two(j) = nf
               end if
            end if
            deallocate(rotated, phi)
         end do
         deallocate(subset, r)
      end do
   end subroutine vss_select_by_keys

   pure function identity_matrix(n) result(a)
      integer, intent(in) :: n !! Order of the square identity matrix.
      real(dp) :: a(n, n)
      integer :: i

      a = 0.0_dp
      do i = 1, n
         a(i, i) = 1.0_dp
      end do
   end function identity_matrix

end module psych_computational_extended
