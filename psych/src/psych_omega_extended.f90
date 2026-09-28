module psych_omega_extended
   use r_kinds, only : dp
   use r_linalg, only : symmetric_eigen
   use psych_correlations, only : squared_multiple_correlations
   use psych_factor, only : procrustes_rotate, promax_rotate
   use psych_factor_extended, only : minimum_residual_factor
   use psych_reliability, only : cronbach_alpha_data
   use psych_reliability_extended, only : omega_statistics
   use psych_remaining, only : ten_berge_bounds
   use psych_stats, only : finite_correlation, finite_mean
   use psych_types, only : alpha_result
   implicit none
   private
   public :: schmid_leiman_decomposition, omega_from_correlation
   public :: direct_schmid_leiman, omega_direct_from_correlation
   public :: glb_factor_reliability, glb_reliability_summary
   public :: test_retest_summary

contains

   subroutine schmid_leiman_decomposition(r, n_factors, loadings, general_loadings, group_loadings, &
      factor_correlations, higher_loadings, communalities, uniqueness, general_fraction, complexity, info, &
      two_factor_option)
      real(dp), intent(in) :: r(:, :) !! Item correlation or covariance matrix to decompose.
      integer, intent(in) :: n_factors !! Number of first-order factors; values are limited to the matrix dimension.
      real(dp), allocatable, intent(out) :: loadings(:, :) !! Orthogonal Schmid-Leiman loadings: general factor then group factors.
      real(dp), allocatable, intent(out) :: general_loadings(:) !! Item loadings on the Schmid-Leiman general factor.
      real(dp), allocatable, intent(out) :: group_loadings(:, :) !! Item loadings on orthogonalized group factors.
      real(dp), allocatable, intent(out) :: factor_correlations(:, :) !! Correlation matrix among oblique first-order factors.
      real(dp), allocatable, intent(out) :: higher_loadings(:) !! First-order factor loadings on the higher-order general factor.
      real(dp), allocatable, intent(out) :: communalities(:) !! First-order common variance for each item.
      real(dp), allocatable, intent(out) :: uniqueness(:) !! Observed item variance minus first-order communality.
      real(dp), allocatable, intent(out) :: general_fraction(:) !! General-factor share of each item's common variance.
      real(dp), allocatable, intent(out) :: complexity(:) !! Hoffmann-style loading complexity for the orthogonal decomposition.
      integer, intent(out) :: info !! Zero on success; nonzero identifies factor extraction or rotation failure.
      integer, intent(in), optional :: two_factor_option !! Two-factor identification: 0 equal, 1 fix first, 2 fix second;
      !! default 0.
      real(dp), allocatable :: orthogonal(:, :), first_uniqueness(:), oblique(:, :), transformation(:, :)
      real(dp), allocatable :: higher(:, :), higher_uniqueness(:)
      real(dp) :: objective, higher_objective, denom
      integer :: i, j, k, local_info, option

      info = 0
      k = min(max(1, n_factors), size(r, 1), size(r, 2))
      option = 0
      if (present(two_factor_option)) option = max(0, min(2, two_factor_option))
      allocate(orthogonal(size(r, 1), k), first_uniqueness(size(r, 1)))
      call minimum_residual_factor(r, k, orthogonal, first_uniqueness, objective, local_info)
      if (local_info /= 0) then
         info = 10 + local_info
         call allocate_empty_schmid(r, k, loadings, general_loadings, group_loadings, factor_correlations, &
            higher_loadings, communalities, uniqueness, general_fraction, complexity)
         return
      end if

      allocate(oblique(size(r, 1), k), transformation(k, k), factor_correlations(k, k))
      if (k == 1) then
         oblique = orthogonal
         transformation = 1.0_dp
         factor_correlations = 1.0_dp
      else
         call promax_rotate(orthogonal, 4, oblique, transformation, factor_correlations, local_info)
         if (local_info /= 0) then
            info = 20 + local_info
            call allocate_empty_schmid(r, k, loadings, general_loadings, group_loadings, factor_correlations, &
               higher_loadings, communalities, uniqueness, general_fraction, complexity)
            return
         end if
      end if

      allocate(higher_loadings(k))
      if (k == 1) then
         higher_loadings = 1.0_dp
      else if (k == 2) then
         select case (option)
         case (1)
            higher_loadings(1) = 1.0_dp
            higher_loadings(2) = factor_correlations(1, 2)
         case (2)
            higher_loadings(1) = factor_correlations(1, 2)
            higher_loadings(2) = 1.0_dp
         case default
            higher_loadings(1) = sqrt(abs(factor_correlations(1, 2)))
            higher_loadings(2) = sign(higher_loadings(1), factor_correlations(1, 2))
         end select
      else
         allocate(higher(k, 1), higher_uniqueness(k))
         call minimum_residual_factor(factor_correlations, 1, higher, higher_uniqueness, higher_objective, local_info)
         if (local_info /= 0) then
            info = 30 + local_info
            call allocate_empty_schmid(r, k, loadings, general_loadings, group_loadings, factor_correlations, &
               higher_loadings, communalities, uniqueness, general_fraction, complexity)
            return
         end if
         higher_loadings = higher(:, 1)
         if (sum(higher_loadings) < 0.0_dp) higher_loadings = -higher_loadings
      end if

      allocate(general_loadings(size(r, 1)), group_loadings(size(r, 1), k), loadings(size(r, 1), k + 1))
      allocate(communalities(size(r, 1)), uniqueness(size(r, 1)), general_fraction(size(r, 1)))
      allocate(complexity(size(r, 1)))
      general_loadings = matmul(oblique, higher_loadings)
      do j = 1, k
         group_loadings(:, j) = oblique(:, j) * sqrt(max(0.0_dp, 1.0_dp - higher_loadings(j) ** 2))
      end do
      loadings(:, 1) = general_loadings
      loadings(:, 2:k + 1) = group_loadings
      communalities = sum(orthogonal * orthogonal, dim=2)
      do i = 1, size(r, 1)
         uniqueness(i) = r(i, i) - communalities(i)
         if (communalities(i) > tiny(1.0_dp)) then
            general_fraction(i) = general_loadings(i) ** 2 / communalities(i)
         else
            general_fraction(i) = 0.0_dp
         end if
         denom = sum(loadings(i, :) ** 4)
         if (denom > tiny(1.0_dp)) then
            complexity(i) = sum(loadings(i, :) ** 2) ** 2 / denom
         else
            complexity(i) = 0.0_dp
         end if
      end do
   end subroutine schmid_leiman_decomposition

   subroutine omega_from_correlation(r, n_factors, omega_h, omega_total, omega_limit, alpha, lambda6, ecv, &
      loadings, group_omega, info, flip_items, two_factor_option)
      real(dp), intent(in) :: r(:, :) !! Item correlation or covariance matrix used for omega estimation.
      integer, intent(in) :: n_factors !! Number of first-order factors in the Schmid-Leiman decomposition.
      real(dp), intent(out) :: omega_h !! Omega hierarchical from the general-factor score variance.
      real(dp), intent(out) :: omega_total !! Omega total from all common-factor loadings.
      real(dp), intent(out) :: omega_limit !! General-factor variance divided by common-score variance.
      real(dp), intent(out) :: alpha !! Cronbach alpha from the keyed item matrix.
      real(dp), intent(out) :: lambda6 !! Guttman lambda-6 from squared multiple correlations.
      real(dp), intent(out) :: ecv !! Explained common variance attributable to the general factor.
      real(dp), allocatable, intent(out) :: loadings(:, :) !! Keyed Schmid-Leiman general and group loadings.
      real(dp), allocatable, intent(out) :: group_omega(:, :) !! Rows general/groups; columns total, general, and group omega.
      integer, intent(out) :: info !! Zero on success; otherwise decomposition status.
      logical, intent(in), optional :: flip_items !! Flip items with negative general loadings; default true as in R `omegah`.
      integer, intent(in), optional :: two_factor_option !! Two-factor identification option passed to Schmid-Leiman.
      real(dp), allocatable :: general(:), groups(:, :), phi(:, :), higher(:), h2(:), u2(:), p2(:), complexity(:)
      real(dp), allocatable :: keyed_r(:, :)
      real(dp) :: common_ss
      integer :: i, k, option
      logical :: do_flip

      option = 0
      if (present(two_factor_option)) option = two_factor_option
      call schmid_leiman_decomposition(r, n_factors, loadings, general, groups, phi, higher, h2, u2, p2, &
         complexity, info, option)
      if (info /= 0) then
         omega_h = 0.0_dp
         omega_total = 0.0_dp
         omega_limit = 0.0_dp
         alpha = 0.0_dp
         lambda6 = 0.0_dp
         ecv = 0.0_dp
         allocate(group_omega(0, 0))
         return
      end if

      allocate(keyed_r(size(r, 1), size(r, 2)))
      keyed_r = r
      do_flip = .true.
      if (present(flip_items)) do_flip = flip_items
      if (do_flip) then
         do i = 1, size(r, 1)
            if (loadings(i, 1) < 0.0_dp) then
               loadings(i, :) = -loadings(i, :)
               keyed_r(i, :) = -keyed_r(i, :)
               keyed_r(:, i) = -keyed_r(:, i)
            end if
         end do
      end if

      call omega_statistics(keyed_r, loadings, omega_h, omega_total, omega_limit, alpha, lambda6)
      common_ss = sum(loadings * loadings)
      if (common_ss > tiny(1.0_dp)) then
         ecv = sum(loadings(:, 1) ** 2) / common_ss
      else
         ecv = 0.0_dp
      end if
      k = size(loadings, 2) - 1
      call omega_group_table(keyed_r, loadings, k, omega_total, group_omega)
   end subroutine omega_from_correlation

   subroutine direct_schmid_leiman(r, n_factors, direct, target, info, cut)
      real(dp), intent(in) :: r(:, :) !! Item correlation matrix for Waller's direct Schmid-Leiman construction.
      integer, intent(in) :: n_factors !! Number of first-order factors used to construct the target.
      real(dp), allocatable, intent(out) :: direct(:, :) !! Direct orthogonal loadings with general factor in column one.
      real(dp), allocatable, intent(out) :: target(:, :) !! General-plus-cluster target used in the Procrustes rotation.
      integer, intent(out) :: info !! Zero when factor extraction, oblique rotation, and Procrustes rotation succeed.
      real(dp), intent(in), optional :: cut !! Minimum absolute oblique loading retained in the cluster target; default 0.3.
      real(dp), allocatable :: unrotated(:, :), uniqueness(:), oblique(:, :), transform(:, :), phi(:, :)
      real(dp), allocatable :: f0(:, :), rotation(:, :)
      real(dp) :: objective, threshold
      integer :: i, j, k, local_info, idx

      k = min(max(1, n_factors), size(r, 1), size(r, 2))
      threshold = 0.3_dp
      if (present(cut)) threshold = max(0.0_dp, cut)
      allocate(unrotated(size(r, 1), k), uniqueness(size(r, 1)))
      call minimum_residual_factor(r, k, unrotated, uniqueness, objective, local_info)
      if (local_info /= 0) then
         info = 10 + local_info
         allocate(direct(0, 0), target(0, 0))
         return
      end if
      allocate(oblique(size(r, 1), k), transform(k, k), phi(k, k))
      if (k == 1) then
         oblique = unrotated
      else
         call promax_rotate(unrotated, 4, oblique, transform, phi, local_info)
         if (local_info /= 0) then
            info = 20 + local_info
            allocate(direct(0, 0), target(0, 0))
            return
         end if
      end if

      allocate(target(size(r, 1), k + 1), f0(size(r, 1), k + 1))
      target = 0.0_dp
      target(:, 1) = 1.0_dp
      do i = 1, size(r, 1)
         idx = maxloc(abs(oblique(i, :)), dim=1)
         if (abs(oblique(i, idx)) > threshold) target(i, idx + 1) = sign(1.0_dp, oblique(i, idx))
      end do
      f0 = 0.0_dp
      f0(:, 2:k + 1) = unrotated
      allocate(direct(size(r, 1), k + 1), rotation(k + 1, k + 1))
      call procrustes_rotate(f0, target, direct, rotation, info)
      if (info /= 0) return
      if (sum(direct(:, 1)) < 0.0_dp) direct(:, 1) = -direct(:, 1)
      do j = 2, k + 1
         if (dot_product(direct(:, j), target(:, j)) < 0.0_dp) direct(:, j) = -direct(:, j)
      end do
   end subroutine direct_schmid_leiman

   subroutine omega_direct_from_correlation(r, n_factors, omega_h, omega_total, omega_limit, alpha, lambda6, ecv, &
      direct, group_omega, info, cut)
      real(dp), intent(in) :: r(:, :) !! Item correlation matrix used by direct Schmid-Leiman omega.
      integer, intent(in) :: n_factors !! Number of first-order group factors.
      real(dp), intent(out) :: omega_h !! General-factor omega hierarchical estimate.
      real(dp), intent(out) :: omega_total !! Total common-factor reliability estimate.
      real(dp), intent(out) :: omega_limit !! General-factor share of reliable variance.
      real(dp), intent(out) :: alpha !! Cronbach alpha from the input matrix.
      real(dp), intent(out) :: lambda6 !! Guttman lambda-6 from the input matrix.
      real(dp), intent(out) :: ecv !! Explained common variance for the direct general factor.
      real(dp), allocatable, intent(out) :: direct(:, :) !! Waller direct Schmid-Leiman loadings.
      real(dp), allocatable, intent(out) :: group_omega(:, :) !! Rows general/groups; columns total, general, and group omega.
      integer, intent(out) :: info !! Zero on successful direct decomposition and reliability calculation.
      real(dp), intent(in), optional :: cut !! Cluster-target loading cutoff; default 0.3.
      real(dp), allocatable :: target(:, :)
      real(dp) :: common_ss, threshold
      integer :: k

      threshold = 0.3_dp
      if (present(cut)) threshold = cut
      call direct_schmid_leiman(r, n_factors, direct, target, info, threshold)
      if (info /= 0) then
         omega_h = 0.0_dp
         omega_total = 0.0_dp
         omega_limit = 0.0_dp
         alpha = 0.0_dp
         lambda6 = 0.0_dp
         ecv = 0.0_dp
         allocate(group_omega(0, 0))
         return
      end if
      call omega_statistics(r, direct, omega_h, omega_total, omega_limit, alpha, lambda6)
      common_ss = sum(direct * direct)
      if (common_ss > tiny(1.0_dp)) then
         ecv = sum(direct(:, 1) ** 2) / common_ss
      else
         ecv = 0.0_dp
      end if
      k = size(direct, 2) - 1
      call omega_group_table(r, direct, k, omega_total, group_omega)
   end subroutine omega_direct_from_correlation

   subroutine glb_factor_reliability(r, glb_value, communalities, n_factors, info)
      real(dp), intent(in) :: r(:, :) !! Item correlation matrix for the factor-based greatest-lower-bound approximation.
      real(dp), intent(out) :: glb_value !! Factor-communality GLB estimate `sum(R with h2 diagonal)/sum(R)`.
      real(dp), allocatable, intent(out) :: communalities(:) !! Communalities from the selected MINRES factor solution.
      integer, intent(out) :: n_factors !! Number of positive reduced-matrix factors retained subject to nonnegative df.
      integer, intent(out) :: info !! Zero on successful SMC, eigensolution, and MINRES extraction.
      real(dp), allocatable :: smc(:), reduced(:, :), values(:), vectors(:, :), loadings(:, :), uniqueness(:), rr(:, :)
      real(dp) :: objective, threshold, total
      integer :: i, p, df, local_info

      p = size(r, 1)
      if (p /= size(r, 2) .or. p < 2) then
         info = 1
         glb_value = 0.0_dp
         n_factors = 0
         allocate(communalities(0))
         return
      end if
      allocate(smc(p), reduced(p, p))
      call squared_multiple_correlations(r, smc, local_info)
      if (local_info /= 0) then
         info = 10 + local_info
         glb_value = 0.0_dp
         n_factors = 0
         allocate(communalities(0))
         return
      end if
      reduced = r
      do i = 1, p
         reduced(i, i) = smc(i)
      end do
      call symmetric_eigen(reduced, values, vectors, local_info, descending=.true.)
      if (local_info /= 0) then
         info = 20 + local_info
         glb_value = 0.0_dp
         n_factors = 0
         allocate(communalities(0))
         return
      end if
      threshold = max(1.0_dp, maxval(abs(values))) * 100.0_dp * epsilon(1.0_dp)
      n_factors = count(values > threshold)
      n_factors = min(max(1, n_factors), p - 1)
      do
         df = p * (p - 1) / 2 - n_factors * p + n_factors * (n_factors - 1) / 2
         if (df >= 0 .or. n_factors <= 1) exit
         n_factors = n_factors - 1
      end do
      allocate(loadings(p, n_factors), uniqueness(p), communalities(p))
      call minimum_residual_factor(r, n_factors, loadings, uniqueness, objective, local_info)
      if (local_info /= 0) then
         info = 30 + local_info
         glb_value = 0.0_dp
         communalities = 0.0_dp
         return
      end if
      communalities = sum(loadings * loadings, dim=2)
      allocate(rr(p, p))
      rr = r
      do i = 1, p
         rr(i, i) = communalities(i)
      end do
      total = sum(r)
      if (abs(total) > tiny(1.0_dp)) then
         glb_value = sum(rr) / total
         info = 0
      else
         glb_value = 0.0_dp
         info = 40
      end if
   end subroutine glb_factor_reliability

   subroutine glb_reliability_summary(r, glb_factor, alpha_pc, r_smc, ten_berge, communalities, n_factors, info)
      real(dp), intent(in) :: r(:, :) !! Item correlation matrix for computational GLB summary measures.
      real(dp), intent(out) :: glb_factor !! Factor-communality GLB approximation corresponding to R `glb.fa`.
      real(dp), intent(out) :: alpha_pc !! Principal-component lower-bound approximation `1 - 1/lambda1`.
      real(dp), intent(out) :: r_smc !! Lambda-6-like reliability based on squared multiple correlations.
      real(dp), intent(out) :: ten_berge(4) !! Ten Berge mu0 through mu3 lower bounds.
      real(dp), allocatable, intent(out) :: communalities(:) !! Communalities used by the factor-based GLB estimate.
      integer, intent(out) :: n_factors !! Number of factors retained for the factor-based GLB estimate.
      integer, intent(out) :: info !! Zero when all eigensolution, SMC, and Ten Berge calculations succeed.
      real(dp), allocatable :: values(:), vectors(:, :), smc(:)
      real(dp) :: total, tr
      integer :: i, local_info

      call glb_factor_reliability(r, glb_factor, communalities, n_factors, info)
      if (info /= 0) then
         alpha_pc = 0.0_dp
         r_smc = 0.0_dp
         ten_berge = 0.0_dp
         return
      end if
      call symmetric_eigen(r, values, vectors, local_info, descending=.true.)
      if (local_info /= 0) then
         info = 50 + local_info
         alpha_pc = 0.0_dp
      else if (size(values) > 0 .and. abs(values(1)) > tiny(1.0_dp)) then
         alpha_pc = 1.0_dp - 1.0_dp / values(1)
      else
         alpha_pc = 0.0_dp
      end if
      allocate(smc(size(r, 1)))
      call squared_multiple_correlations(r, smc, local_info)
      if (local_info /= 0) then
         info = max(info, 60 + local_info)
         r_smc = 0.0_dp
      else
         total = sum(r)
         tr = 0.0_dp
         do i = 1, size(r, 1)
            tr = tr + r(i, i)
         end do
         if (abs(total) > tiny(1.0_dp)) then
            r_smc = (total + sum(smc) - tr) / total
         else
            r_smc = 0.0_dp
         end if
      end if
      call ten_berge_bounds(r, ten_berge(1), ten_berge(2), ten_berge(3), ten_berge(4), local_info)
      if (local_info /= 0) info = max(info, 70 + local_info)
   end subroutine glb_reliability_summary

   subroutine test_retest_summary(time1, time2, scale_correlation, raw_alpha, standardized_alpha, g6, &
      item_correlations, profile_correlations, info)
      real(dp), intent(in) :: time1(:, :) !! Subject-by-item measurements at the first occasion; NaNs are allowed.
      real(dp), intent(in) :: time2(:, :) !! Matched subject-by-item measurements at the second occasion; NaNs are allowed.
      real(dp), intent(out) :: scale_correlation !! Correlation of subject mean scale scores across the two occasions.
      real(dp), intent(out) :: raw_alpha(2) !! Raw Cronbach alpha for time one and time two.
      real(dp), intent(out) :: standardized_alpha(2) !! Standardized Cronbach alpha for time one and time two.
      real(dp), intent(out) :: g6(2) !! Guttman lambda-6 for time one and time two.
      real(dp), allocatable, intent(out) :: item_correlations(:) !! Matched-time correlation for each item over subjects.
      real(dp), allocatable, intent(out) :: profile_correlations(:) !! Within-subject correlation of item profiles across time.
      integer, intent(out) :: info !! Zero on compatible dimensions and successful alpha computations.
      real(dp), allocatable :: score1(:), score2(:)
      type(alpha_result) :: a1, a2
      integer :: i, p, n, info1, info2

      if (any(shape(time1) /= shape(time2))) then
         scale_correlation = 0.0_dp
         raw_alpha = 0.0_dp
         standardized_alpha = 0.0_dp
         g6 = 0.0_dp
         allocate(item_correlations(0), profile_correlations(0))
         info = 1
         return
      end if
      n = size(time1, 1)
      p = size(time1, 2)
      allocate(score1(n), score2(n), item_correlations(p), profile_correlations(n))
      do i = 1, n
         score1(i) = finite_mean(time1(i, :))
         score2(i) = finite_mean(time2(i, :))
         profile_correlations(i) = finite_correlation(time1(i, :), time2(i, :))
      end do
      scale_correlation = finite_correlation(score1, score2)
      do i = 1, p
         item_correlations(i) = finite_correlation(time1(:, i), time2(:, i))
      end do
      call cronbach_alpha_data(time1, a1, info1)
      call cronbach_alpha_data(time2, a2, info2)
      raw_alpha = [a1%raw_alpha, a2%raw_alpha]
      standardized_alpha = [a1%standardized_alpha, a2%standardized_alpha]
      g6 = [a1%g6, a2%g6]
      info = max(info1, info2)
   end subroutine test_retest_summary

   subroutine omega_group_table(r, loadings, n_groups, omega_total, table)
      real(dp), intent(in) :: r(:, :) !! Keyed item correlation matrix.
      real(dp), intent(in) :: loadings(:, :) !! General factor in column one followed by group-factor loadings.
      integer, intent(in) :: n_groups !! Number of group-factor columns after the general factor.
      real(dp), intent(in) :: omega_total !! Overall omega total placed in row one of the table.
      real(dp), allocatable, intent(out) :: table(:, :) !! Rows general/groups; columns total, general, and group omega.
      logical, allocatable :: selected(:)
      real(dp) :: general_sum, group_sum, total_variance, v_group, group_sum_sq
      integer :: g, i, j, idx

      allocate(table(n_groups + 1, 3), selected(size(r, 1)))
      table = 0.0_dp
      total_variance = sum(r)
      group_sum_sq = 0.0_dp
      do g = 1, n_groups
         selected = .false.
         do i = 1, size(r, 1)
            idx = maxloc(abs(loadings(i, 2:n_groups + 1)), dim=1)
            if (idx == g) selected(i) = .true.
         end do
         if (.not. any(selected)) cycle
         v_group = 0.0_dp
         do j = 1, size(r, 2)
            if (.not. selected(j)) cycle
            do i = 1, size(r, 1)
               if (selected(i)) v_group = v_group + r(i, j)
            end do
         end do
         group_sum = sum(loadings(:, g + 1), mask=selected)
         general_sum = sum(loadings(:, 1), mask=selected)
         group_sum_sq = group_sum_sq + group_sum ** 2
         if (v_group > tiny(1.0_dp)) then
            table(g + 1, 1) = (group_sum ** 2 + general_sum ** 2) / v_group
            table(g + 1, 2) = general_sum ** 2 / v_group
            table(g + 1, 3) = group_sum ** 2 / v_group
         end if
      end do
      table(1, 1) = omega_total
      if (abs(total_variance) > tiny(1.0_dp)) then
         table(1, 2) = sum(loadings(:, 1)) ** 2 / total_variance
         table(1, 3) = group_sum_sq / total_variance
      end if
   end subroutine omega_group_table

   subroutine allocate_empty_schmid(r, k, loadings, general, groups, phi, higher, h2, u2, p2, complexity)
      real(dp), intent(in) :: r(:, :) !! Matrix whose row count determines empty-output item dimension.
      integer, intent(in) :: k !! Requested number of first-order factors.
      real(dp), allocatable, intent(out) :: loadings(:, :) !! Zeroed general-plus-group loading matrix.
      real(dp), allocatable, intent(out) :: general(:) !! Zeroed general loadings.
      real(dp), allocatable, intent(out) :: groups(:, :) !! Zeroed group loadings.
      real(dp), allocatable, intent(out) :: phi(:, :) !! Zeroed factor-correlation matrix.
      real(dp), allocatable, intent(out) :: higher(:) !! Zeroed higher-order factor loadings.
      real(dp), allocatable, intent(out) :: h2(:) !! Zeroed communalities.
      real(dp), allocatable, intent(out) :: u2(:) !! Zeroed uniquenesses.
      real(dp), allocatable, intent(out) :: p2(:) !! Zeroed general common-variance fractions.
      real(dp), allocatable, intent(out) :: complexity(:) !! Zeroed item complexities.

      if (allocated(loadings)) deallocate(loadings)
      if (allocated(general)) deallocate(general)
      if (allocated(groups)) deallocate(groups)
      if (allocated(phi)) deallocate(phi)
      if (allocated(higher)) deallocate(higher)
      if (allocated(h2)) deallocate(h2)
      if (allocated(u2)) deallocate(u2)
      if (allocated(p2)) deallocate(p2)
      if (allocated(complexity)) deallocate(complexity)
      allocate(loadings(size(r, 1), k + 1), general(size(r, 1)), groups(size(r, 1), k))
      allocate(phi(k, k), higher(k), h2(size(r, 1)), u2(size(r, 1)), p2(size(r, 1)), complexity(size(r, 1)))
      loadings = 0.0_dp
      general = 0.0_dp
      groups = 0.0_dp
      phi = 0.0_dp
      higher = 0.0_dp
      h2 = 0.0_dp
      u2 = 0.0_dp
      p2 = 0.0_dp
      complexity = 0.0_dp
   end subroutine allocate_empty_schmid

end module psych_omega_extended
