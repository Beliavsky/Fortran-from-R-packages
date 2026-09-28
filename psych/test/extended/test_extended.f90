program test_extended
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
   use psych
   implicit none

   integer :: failures

   failures = 0
   call test_inference(failures)
   call test_factor_extensions(failures)
   call test_scoring_and_irt(failures)
   call test_validity_and_vss(failures)
   call test_additional_translations(failures)
   call test_simulation_algorithms(failures)
   call test_omega_extensions(failures)
   call test_parallel_extensions(failures)
   call test_selection_extensions(failures)
   call test_iclust_extensions(failures)

   if (failures /= 0) then
      error stop "psych extended tests failed"
   end if
   print '(a)', "All psych extended tests passed"

contains

   subroutine assert_close(actual, expected, tolerance, message, failures)
      real(dp), intent(in) :: actual !! Actual scalar value produced by the routine under test.
      real(dp), intent(in) :: expected !! Reference scalar value.
      real(dp), intent(in) :: tolerance !! Maximum accepted absolute error.
      character(len=*), intent(in) :: message !! Failure description printed when the assertion fails.
      integer, intent(inout) :: failures !! Running number of failed assertions.

      if (abs(actual - expected) > tolerance) then
         print '(a,2(1x,es14.6))', trim(message), actual, expected
         failures = failures + 1
      end if
   end subroutine assert_close

   subroutine assert_true(condition, message, failures)
      logical, intent(in) :: condition !! Boolean condition that must be true.
      character(len=*), intent(in) :: message !! Failure description printed when the assertion fails.
      integer, intent(inout) :: failures !! Running number of failed assertions.

      if (.not. condition) then
         print '(a)', trim(message)
         failures = failures + 1
      end if
   end subroutine assert_true

   subroutine test_inference(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      real(dp) :: lower, upper, rbar, statistic, p_value, df, chisq, rmsea, srmr
      real(dp) :: r(3, 3)
      integer :: info

      call alpha_confidence_interval(0.8_dp, 100, 10, 0.05_dp, lower, upper, rbar)
      call assert_true(lower < 0.8_dp .and. upper > 0.8_dp, "alpha confidence interval brackets estimate", failures)
      call assert_true(rbar > 0.0_dp .and. rbar < 1.0_dp, "alpha confidence interval returns mean r", failures)

      call effect_size_confidence_interval(0.5_dp, 50, 50, 0.05_dp, lower, upper)
      call assert_true(lower < 0.5_dp .and. upper > 0.5_dp, "effect-size interval brackets d", failures)

      call correlation_test_general(100, 0.3_dp, statistic, p_value, df, r34=0.3_dp, n2=100, pooled=.false.)
      call assert_close(statistic, 0.0_dp, 1.0e-12_dp, "independent equal correlations z", failures)
      call assert_close(p_value, 1.0_dp, 1.0e-10_dp, "independent equal correlations p", failures)

      r = reshape([1.0_dp, 0.3_dp, 0.2_dp, 0.3_dp, 1.0_dp, 0.25_dp, 0.2_dp, 0.25_dp, 1.0_dp], [3, 3])
      call jennrich_correlation_test(r, r, 100, 120, chisq, p_value, df, info)
      call assert_true(info == 0, "Jennrich inversion succeeds", failures)
      call assert_close(chisq, 0.0_dp, 1.0e-10_dp, "Jennrich identical matrices", failures)
      call correlation_normal_test(r, 100, .true., chisq, p_value, df, rmsea, srmr, r2=r, n2=100)
      call assert_close(chisq, 0.0_dp, 1.0e-10_dp, "normal correlation test identical matrices", failures)
   end subroutine test_inference

   subroutine test_factor_extensions(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      real(dp) :: r(4, 4), objective, rms, fit, chisq, df, p_value
      real(dp), allocatable :: loadings(:, :), uniqueness(:), beta(:, :), residual(:, :), r2(:)
      real(dp), allocatable :: rotated(:, :), rotation(:, :), target(:, :), discrimination(:, :), difficulty(:, :, :)
      real(dp) :: thresholds(4, 1)
      integer :: info

      r = 0.5_dp
      r = r + 0.5_dp * reshape([1.0_dp,0.0_dp,0.0_dp,0.0_dp, &
         0.0_dp,1.0_dp,0.0_dp,0.0_dp, 0.0_dp,0.0_dp,1.0_dp,0.0_dp, &
         0.0_dp,0.0_dp,0.0_dp,1.0_dp], [4,4])
      call minimum_residual_factor(r, 1, loadings, uniqueness, objective, info)
      call assert_true(info == 0, "MINRES eigensolutions succeed", failures)
      call assert_true(objective < 0.05_dp, "MINRES equicorrelation residual is small", failures)
      call factor_fit_statistics(r, loadings, 200, rms, fit, chisq, df, p_value)
      call assert_true(rms < 0.15_dp .and. fit > 0.8_dp, "factor fit statistics reflect good one-factor fit", failures)

      call factor_regression_partition(r, [1,2], [3,4], beta, residual, r2, info)
      call assert_true(info == 0, "factor regression partition inversion succeeds", failures)
      call assert_true(all(r2 > 0.0_dp), "factor regression returns positive multiple R2", failures)

      call equamax_rotate(reshape([0.8_dp,0.7_dp,0.2_dp,0.1_dp,0.1_dp,0.2_dp,0.7_dp,0.8_dp], [4,2]), &
         rotated, rotation, info)
      call assert_true(info == 0, "equamax rotation succeeds", failures)
      call assert_close(sum(rotated * rotated), 2.36_dp, 1.0e-8_dp, "equamax preserves total communality", failures)

      allocate(target(4,2))
      target = reshape([0.8_dp,0.7_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.7_dp,0.8_dp], [4,2])
      call target_orthogonal_rotate(rotated, target, loadings, rotation, info)
      call assert_true(info == 0, "orthogonal target rotation succeeds", failures)

      thresholds(:,1) = 0.0_dp
      call factor_to_irt_parameters(reshape([0.5_dp,0.4_dp,0.3_dp,0.2_dp], [4,1]), thresholds, &
         discrimination, difficulty)
      call assert_close(discrimination(1,1), 0.5_dp / sqrt(0.75_dp), 1.0e-12_dp, &
         "factor-to-IRT discrimination", failures)
   end subroutine test_factor_extensions

   subroutine test_scoring_and_irt(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      integer :: clusters(4), items(4,2), levels(2)
      real(dp) :: x(4,2), theta_in(4), loadings(2,1), difficulty_in(2,1), discrimination_in(2)
      real(dp), allocatable :: keys(:, :), means(:, :), sds(:, :), thresholds(:, :), bin_centers(:)
      real(dp), allocatable :: proportions(:, :, :), theta(:), fit(:)
      integer, allocatable :: counts(:, :), group_values(:)

      clusters = [1,1,2,2]
      call cluster_assignments_to_keys(clusters, keys)
      call assert_true(size(keys,2) == 2 .and. abs(sum(keys) - 4.0_dp) < 1.0e-12_dp, &
         "cluster2keys translation", failures)

      x = reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp, 2.0_dp,4.0_dp,6.0_dp,8.0_dp], [4,2])
      call grouped_describe_matrix(x, [1,1,2,2], group_values, means, sds, counts)
      call assert_close(means(1,1), 1.5_dp, 1.0e-12_dp, "grouped describe first group", failures)
      call assert_close(means(2,2), 7.0_dp, 1.0e-12_dp, "grouped describe second group", failures)

      items = reshape([0,0,1,1, 0,1,0,1], [4,2])
      call irt_thresholds(items, 0, 1, thresholds)
      call assert_close(thresholds(1,1), 0.0_dp, 1.0e-8_dp, "balanced dichotomous threshold", failures)
      theta_in = [-1.5_dp,-0.5_dp,0.5_dp,1.5_dp]
      levels = [0,1]
      call irt_response_curves(theta_in, items, levels, 2, bin_centers, proportions, counts)
      call assert_close(sum(proportions(1,1,:)), 1.0_dp, 1.0e-12_dp, "response-curve proportions sum to one", failures)

      loadings(:,1) = [0.5_dp,0.6_dp]
      difficulty_in(:,1) = 0.0_dp
      discrimination_in = [1.0_dp,1.0_dp]
      call irt_score_poly(items, 0, difficulty_in, discrimination_in, -4.0_dp, 4.0_dp, theta, fit)
      call assert_true(theta(4) > theta(1), "ordinal IRT score increases with positive responses", failures)
   end subroutine test_scoring_and_irt

   subroutine test_validity_and_vss(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      real(dp) :: items(6,3), criteria(6,1), keys(3,2), q(2,2), margins(2), r(4,4)
      real(dp) :: xmat(2,2), ymat(1,1), xymat(2,1), set_r2, u, alpha_std, average_r, median_r, map_value, f1_f2
      real(dp), allocatable :: validity(:, :), predicted(:, :), asymptotic(:, :), tetra(:, :), phi(:, :)
      real(dp), allocatable :: beta(:, :), r2(:), residual(:, :), fit(:), rss(:)
      real(dp) :: loadings(4,2), factor_phi(2,2), rho, se, lower, upper

      items(:,1) = [1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp]
      items(:,2) = items(:,1)
      items(:,3) = -items(:,1)
      criteria(:,1) = items(:,1)
      keys = reshape([1.0_dp,1.0_dp,0.0_dp, 1.0_dp,0.0_dp,-1.0_dp], [3,2])
      call item_validity_by_keys(items, criteria, keys, validity)
      call assert_close(validity(1,1), 1.0_dp, 1.0e-12_dp, "item validity first key", failures)
      call assert_close(validity(2,1), 1.0_dp, 1.0e-12_dp, "item validity reversed key", failures)
      call predicted_validity_summary(validity, [0.5_dp,0.5_dp], [2,2], predicted, asymptotic)
      call assert_true(all(predicted > 0.0_dp) .and. all(asymptotic > predicted), "predicted validity scaling", failures)

      call yule_bonett_interval(reshape([40.0_dp,10.0_dp,10.0_dp,40.0_dp],[2,2]), 1.0_dp, .false., &
         0.05_dp, rho, se, lower, upper)
      call assert_true(rho > 0.0_dp .and. lower < rho .and. upper > rho, "Yule-Bonett interval", failures)
      q = reshape([1.0_dp,0.5_dp,0.5_dp,1.0_dp], [2,2])
      margins = [0.5_dp,0.5_dp]
      call yule_to_phi_values(q, margins, phi)
      call yule_to_tetrachoric_values(q, margins, tetra)
      call assert_true(phi(1,2) > 0.0_dp .and. tetra(1,2) > phi(1,2), "Yule transforms positive", failures)

      xmat = reshape([1.0_dp,0.2_dp,0.2_dp,1.0_dp], [2,2])
      ymat(1,1) = 1.0_dp
      xymat(:,1) = [0.5_dp,0.3_dp]
      call regression_from_matrices(xmat, ymat, xymat, beta, r2, residual, set_r2)
      call assert_true(r2(1) > 0.25_dp .and. r2(1) < 0.35_dp, "matrix regression R2", failures)
      call assert_true(set_r2 > 0.0_dp .and. set_r2 < 1.0_dp, "set correlation determinant reduction", failures)

      r = 0.5_dp
      r = r + 0.5_dp * reshape([1.0_dp,0.0_dp,0.0_dp,0.0_dp, &
         0.0_dp,1.0_dp,0.0_dp,0.0_dp, 0.0_dp,0.0_dp,1.0_dp,0.0_dp, &
         0.0_dp,0.0_dp,0.0_dp,1.0_dp], [4,4])
      call unidimensionality_statistics(r, 2, u, alpha_std, average_r, median_r, map_value, f1_f2)
      call assert_true(u > 0.7_dp .and. alpha_std > 0.7_dp, "unidimensionality equicorrelation", failures)
      loadings = reshape([0.7_dp,0.7_dp,0.1_dp,0.1_dp, 0.1_dp,0.1_dp,0.7_dp,0.7_dp], [4,2])
      factor_phi = 0.0_dp
      factor_phi(1,1) = 1.0_dp
      factor_phi(2,2) = 1.0_dp
      call vss_complexity_fit(r, loadings, factor_phi, 2, fit, rss)
      call assert_true(fit(2) >= fit(1), "VSS fit does not degrade at greater complexity", failures)
   end subroutine test_validity_and_vss


   subroutine test_additional_translations(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      real(dp) :: prep, dprime, probability, requiv
      real(dp) :: xprofile(5, 2), fload(2, 1), rgrice(2, 2)
      real(dp) :: original(4, 4), loadings(4, 1), clusters(4, 1)
      real(dp) :: cluster_fit, structure_fit, pattern_fit, cluster_rmse, pattern_rmse
      real(dp) :: xblock(2, 2), xyblock(1, 2), yblock(1, 1)
      real(dp) :: xdata(4, 2), ydata(4, 1), r12(2, 2)
      real(dp), allocatable :: similarity(:, :), r2(:), h(:, :), combined(:, :), cor2(:, :)
      real(dp), allocatable :: a1(:, :), a2(:, :), scale_residual(:, :), scale_probability(:, :)
      real(dp), allocatable :: congeneric(:, :), pattern(:, :), hierarchical(:, :), dichot(:, :)
      real(dp), allocatable :: best_period(:), cv(:), cv_matrix(:, :)
      real(dp), allocatable :: eigenvalues(:), normal_lower(:), normal_upper(:), boot_lower(:), boot_upper(:)
      real(dp), allocatable :: count_quantiles(:)
      real(dp), allocatable :: thurstone_scale(:), thurstone_choice(:, :), thurstone_residual(:, :)
      real(dp), allocatable :: polyserial(:, :), general_model(:, :), bonds_cor(:, :), bonds_loadings(:, :)
      real(dp), allocatable :: vss_items(:, :)
      real(dp), allocatable :: npl_items(:, :), npl_theta(:), npn_items(:, :), npn_theta(:), poly_marginal(:, :)
      real(dp) :: gof, original_ss, residual_ss, com_table(2, 2), phi_value, yule_value, tetra_value
      real(dp) :: grm_sum, angle_period(48), yperiod(48, 1), x_eigen(20, 2), x_complete(6, 4)
      real(dp) :: cv_weights(1, 1), cv_cor(2, 2)
      real(dp) :: ten_r(3, 3), mu0, mu1, mu2, mu3, thurstone_gof
      real(dp) :: thurstone_input(2, 2), poly_x(4, 1), poly_y(4, 1)
      real(dp) :: latent_cor2(2, 2), marginal_thresholds(2, 2)
      type(cosinor_result), allocatable :: best_cos(:)
      type(describe_result) :: pair_desc
      integer, allocatable :: design(:, :), contingency(:, :)
      integer, allocatable :: valid_counts(:), detected_types(:), group_values(:), group_counts(:, :)
      integer :: info, i, j, duplicates, complete_cases
      integer :: mixed_types(3), group_labels(6)
      real(dp) :: data_summary(4, 3), mixed_data(10, 3), group_data(6, 1)
      real(dp), allocatable :: head_values(:, :), tail_values(:, :), big_r(:, :), mixed_r(:, :), polydi_r(:, :)
      real(dp), allocatable :: group_means(:, :), group_sds(:, :), grand_means(:), ms_between(:), ms_within(:)
      real(dp), allocatable :: f_values(:), icc1(:), icc2(:)
      logical :: pattern_available

      call replication_probability(0.05_dp, prep, n=100, d_prime=dprime, r_equiv=requiv)
      call assert_close(prep, 0.8776029281733251_dp, 2.0e-7_dp, "p.rep probability", failures)
      call assert_close(dprime, 0.4009229697505582_dp, 2.0e-6_dp, "p.rep d-prime", failures)
      call assert_close(requiv, 0.1965511955932256_dp, 2.0e-6_dp, "p.rep equivalent r", failures)

      call replication_probability_f(4.0_dp, 98.0_dp, prep, dprime, probability, requiv)
      call assert_close(prep, 0.8800288807502004_dp, 2.0e-7_dp, "p.rep.f probability", failures)
      call assert_close(probability, 0.0482677719493676_dp, 2.0e-7_dp, "p.rep.f tail probability", failures)

      call replication_probability_r(0.3_dp, 100, prep, dprime, probability)
      call assert_close(prep, 0.9853274459125791_dp, 2.0e-7_dp, "p.rep.r probability", failures)
      call assert_close(probability, 0.0020608812020897_dp, 2.0e-7_dp, "p.rep.r two-tailed probability", failures)

      call replication_probability_t(2.0_dp, 30.0_dp, prep, dprime, probability, requiv)
      call assert_close(prep, 0.8712842968137269_dp, 2.0e-7_dp, "p.rep.t probability", failures)
      call assert_close(dprime, 0.7302967433402214_dp, 2.0e-12_dp, "p.rep.t d-prime", failures)

      xprofile(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
      xprofile(:, 2) = [5.0_dp, 4.0_dp, 3.0_dp, 2.0_dp, 1.0_dp]
      call cohen_profile_similarity(xprofile, similarity)
      call assert_close(similarity(1, 1), 1.0_dp, 1.0e-12_dp, "Cohen profile self similarity", failures)
      call assert_close(similarity(2, 1), -1.0_dp, 1.0e-12_dp, "Cohen profile reverse similarity", failures)

      fload(:, 1) = [0.8_dp, 0.6_dp]
      call factor_score_indeterminacy(fload, r2, info, h=h)
      call assert_true(info == 0, "factor score indeterminacy succeeds", failures)
      call assert_close(r2(1), 0.7006237006237006_dp, 2.0e-8_dp, "factor score indeterminacy R2", failures)
      rgrice = reshape([1.0_dp, 0.3_dp, 0.3_dp, 1.0_dp], [2, 2])
      call factor_score_indeterminacy(fload, r2, info, r=rgrice, grice=.true.)
      call assert_close(r2(1), 0.7824175824175824_dp, 2.0e-8_dp, "Grice factor score indeterminacy", failures)

      original = 0.4_dp
      do i = 1, 4
         original(i, i) = 1.0_dp
      end do
      loadings(:, 1) = sqrt(0.4_dp)
      clusters(:, 1) = 1.0_dp
      call cluster_fit_metrics(original, loadings, clusters, cluster_fit, structure_fit, pattern_fit, &
         cluster_rmse, pattern_rmse, pattern_available, info)
      call assert_true(info == 0 .and. pattern_available, "cluster.fit pattern model available", failures)
      call assert_close(cluster_fit, 0.7567567567567568_dp, 2.0e-10_dp, "cluster.fit cluster fit", failures)
      call assert_close(structure_fit, 0.7567567567567568_dp, 2.0e-10_dp, "cluster.fit structure fit", failures)
      call assert_close(pattern_fit, 1.0_dp, 2.0e-10_dp, "cluster.fit pattern fit", failures)
      call assert_close(cluster_rmse, 0.3464101615137755_dp, 2.0e-10_dp, "cluster.fit RMSE", failures)
      call assert_close(pattern_rmse, 0.0_dp, 1.0e-8_dp, "cluster.fit pattern RMSE", failures)

      xblock = reshape([1.0_dp, 0.2_dp, 0.2_dp, 1.0_dp], [2, 2])
      xyblock(1, :) = [0.4_dp, 0.3_dp]
      yblock(1, 1) = 1.0_dp
      call super_correlation_matrix(xblock, xyblock, combined, y=yblock, info=info)
      call assert_true(info == 0 .and. all(shape(combined) == [3, 3]), "superCor block dimensions", failures)
      call assert_close(combined(3, 1), 0.4_dp, 1.0e-12_dp, "superCor lower-left block", failures)
      call assert_close(combined(1, 3), 0.4_dp, 1.0e-12_dp, "superCor upper-right block", failures)

      xdata(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
      xdata(:, 2) = [2.0_dp, 4.0_dp, 6.0_dp, 8.0_dp]
      ydata(:, 1) = [4.0_dp, 3.0_dp, 2.0_dp, 1.0_dp]
      call cor2_matrix(xdata, cor2, y=ydata)
      call assert_close(cor2(1, 1), -1.0_dp, 1.0e-12_dp, "cor2 cross correlation", failures)
      call assert_close(cor2(2, 1), -1.0_dp, 1.0e-12_dp, "cor2 second cross correlation", failures)

      r12 = reshape([0.8_dp, 0.0_dp, 0.0_dp, 0.4_dp], [2, 2])
      call interbattery_loadings(r12, 2, a1, a2, info)
      call assert_true(info == 0, "interbattery eigensolution", failures)
      call assert_close(a1(1, 1), sqrt(0.8_dp), 2.0e-10_dp, "interbattery first X loading", failures)
      call assert_close(a2(1, 1), sqrt(0.8_dp), 2.0e-10_dp, "interbattery first Y loading", failures)
      call assert_close(a1(2, 2), sqrt(0.4_dp), 2.0e-10_dp, "interbattery second X loading", failures)



      grm_sum = 0.0_dp
      do i = 1, 5
         grm_sum = grm_sum + graded_response_probability(0.0_dp, i, [-1.5_dp, -0.5_dp, 0.5_dp, 1.5_dp])
      end do
      call assert_close(grm_sum, 1.0_dp, 1.0e-12_dp, "logistic.grm category probabilities sum to one", failures)
      call assert_close(graded_response_probability(0.0_dp, 3, [-1.5_dp, -0.5_dp, 0.5_dp, 1.5_dp]), &
         0.3583573983507859_dp, 2.0e-12_dp, "logistic.grm middle category", failures)

      call integer_contingency_table([1, 1, 2, 2, 2], [0, 1, 0, 1, 1], contingency)
      call assert_true(all(shape(contingency) == [2, 2]), "tableF dimensions", failures)
      call assert_true(all(contingency == reshape([1, 1, 1, 2], [2, 2])), "tableF counts", failures)

      call scaling_fit([0.0_dp, 1.0_dp], reshape([0.5_dp, 0.2689414213699951_dp, &
         0.7310585786300049_dp, 0.5_dp], [2, 2]), gof, original_ss, residual_ss, scale_residual, &
         probability=scale_probability)
      call assert_close(gof, 1.0_dp, 2.0e-12_dp, "scaling.fits exact logit fit", failures)
      call assert_close(residual_ss, 0.0_dp, 2.0e-24_dp, "scaling.fits residual sum of squares", failures)
      call assert_close(scale_probability(1, 2), 0.7310585786300049_dp, 2.0e-12_dp, &
         "scaling.fits rowwise probability", failures)

      call make_congeneric_model(model=congeneric, pattern=pattern)
      call assert_true(all(shape(congeneric) == [4, 4]) .and. all(shape(pattern) == [4, 5]), &
         "make.congeneric dimensions", failures)
      call assert_close(congeneric(1, 1), 1.0_dp, 1.0e-12_dp, "make.congeneric unit diagonal", failures)
      call assert_close(congeneric(1, 2), 0.56_dp, 1.0e-12_dp, "make.congeneric off diagonal", failures)

      call make_hierarchical_model(hierarchical, info=info)
      call assert_true(info == 0 .and. all(shape(hierarchical) == [9, 9]), "make.hierarchical dimensions", failures)
      call assert_close(hierarchical(1, 2), 0.56_dp, 1.0e-12_dp, "make.hierarchical within factor", failures)
      call assert_close(hierarchical(1, 4), 0.4032_dp, 1.0e-12_dp, "make.hierarchical cross factor", failures)

      call simulate_item_dichot(8, 6, dichot, xloading=0.0_dp, yloading=0.0_dp, gloading=1.0_dp, &
         low=-huge(1.0_dp), high=huge(1.0_dp))
      call assert_true(all(shape(dichot) == [6, 8]), "item.dichot dimensions", failures)
      call assert_true(minval(dichot(:, 1:4)) > 0.5_dp .and. maxval(dichot(:, 5:8)) < 0.5_dp, &
         "item.dichot threshold regrouping", failures)

      call comorbidity_statistics(0.2_dp, 0.3_dp, 0.1_dp, com_table, phi_value, yule_value, tetra_value)
      call assert_close(sum(com_table), 1.0_dp, 1.0e-12_dp, "comorbidity probability table", failures)
      call assert_true(phi_value > 0.0_dp .and. yule_value > 0.0_dp .and. tetra_value > 0.0_dp, &
         "comorbidity association measures", failures)

      do i = 1, 48
         angle_period(i) = real(i - 1, dp)
         yperiod(i, 1) = cos(2.0_dp * acos(-1.0_dp) * angle_period(i) / 24.0_dp)
      end do
      call cosinor_period_fit(angle_period, yperiod, [22.0_dp, 23.0_dp, 24.0_dp, 25.0_dp, 26.0_dp], &
         best_period, best_cos, hours=.true.)
      call assert_close(best_period(1), 24.0_dp, 1.0e-12_dp, "cosinor.period selects generating period", failures)
      call assert_true(best_cos(1)%fit > 0.999999_dp, "cosinor.period exact cosine fit", failures)

      cv_weights(1, 1) = 1.0_dp
      cv_cor = reshape([1.0_dp, 0.5_dp, 0.5_dp, 1.0_dp], [2, 2])
      call cross_validation_from_correlation(cv_weights, cv_cor, cv, cv_matrix, info)
      call assert_true(info == 0, "crossValidation correlation dimensions", failures)
      call assert_close(cv(1), 0.5_dp, 1.0e-12_dp, "crossValidation criterion correlation", failures)
      call assert_close(cv_matrix(1, 1), 0.5_dp, 1.0e-12_dp, "crossValidation matrix correlation", failures)

      do i = 1, 20
         x_eigen(i, 1) = real(i, dp)
         x_eigen(i, 2) = real(i * i, dp)
      end do
      call eigen_confidence_intervals(x_eigen, eigenvalues, normal_lower, normal_upper, boot_lower, boot_upper, &
         n_iter=20, info=info)
      call assert_true(info == 0 .and. size(eigenvalues) == 2, "eigenCi dimensions", failures)
      call assert_true(all(normal_lower <= eigenvalues) .and. all(eigenvalues <= normal_upper), &
         "eigenCi normal intervals bracket eigenvalues", failures)
      call assert_true(all(boot_lower <= boot_upper), "eigenCi bootstrap interval ordering", failures)

      do j = 1, 4
         do i = 1, 6
            x_complete(i, j) = real(i + j, dp)
         end do
      end do
      call pairwise_sample_description(x_complete, pair_desc, count_quantiles, sample_size=3)
      call assert_close(pair_desc%mean, 6.0_dp, 1.0e-12_dp, "pairwiseSample complete-data mean count", failures)
      call assert_true(all(abs(count_quantiles - 6.0_dp) < 1.0e-12_dp), &
         "pairwiseSample complete-data quantiles", failures)

      ten_r = 0.5_dp
      do i = 1, 3
         ten_r(i, i) = 1.0_dp
      end do
      call ten_berge_bounds(ten_r, mu0, mu1, mu2, mu3, info)
      call assert_true(info == 0, "tenberge succeeds", failures)
      call assert_close(mu0, 0.75_dp, 1.0e-12_dp, "tenberge mu0", failures)
      call assert_close(mu1, 0.75_dp, 1.0e-12_dp, "tenberge mu1", failures)
      call assert_close(mu2, 0.75_dp, 1.0e-12_dp, "tenberge mu2", failures)
      call assert_close(mu3, 0.75_dp, 1.0e-12_dp, "tenberge mu3", failures)

      thurstone_input = reshape([0.5_dp, 0.25_dp, 0.75_dp, 0.5_dp], [2, 2])
      call thurstone_scaling(thurstone_input, thurstone_scale, thurstone_gof, thurstone_choice, &
         thurstone_residual, info=info)
      call assert_true(info == 0, "thurstone succeeds", failures)
      call assert_close(thurstone_scale(1), 0.0_dp, 1.0e-10_dp, "thurstone first scale", failures)
      call assert_close(thurstone_scale(2), 0.6744897501960817_dp, 2.0e-7_dp, "thurstone second scale", failures)
      call assert_close(thurstone_gof, 1.0_dp, 2.0e-7_dp, "thurstone exact fit", failures)
      call assert_true(maxval(abs(thurstone_residual)) < 2.0e-7_dp, "thurstone residual", failures)

      poly_x(:, 1) = [-1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
      poly_y(:, 1) = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
      call polyserial_correlations(poly_x, poly_y, polyserial, info)
      call assert_true(info == 0, "polyserial succeeds", failures)
      call assert_close(polyserial(1, 1), 0.8862269254527580_dp, 2.0e-7_dp, "polyserial binary approximation", failures)

      call simulate_general_structure(6, 2, 0.2_dp, 0.3_dp, general_model, info=info)
      call assert_true(info == 0 .and. all(shape(general_model) == [6, 6]), "sim.general dimensions", failures)
      call assert_close(general_model(1, 2), 0.5_dp, 1.0e-12_dp, "sim.general within group", failures)
      call assert_close(general_model(1, 4), 0.2_dp, 1.0e-12_dp, "sim.general across groups", failures)

      call simulate_bonds_model(5, [0.0_dp, 0.0_dp, 0.5_dp, 0.6_dp], 0.8_dp, bonds_cor, bonds_loadings, info=info)
      call assert_true(info == 0 .and. all(shape(bonds_loadings) == [9, 4]), "sim.bonds loading dimensions", failures)
      call assert_true(all(shape(bonds_cor) == [9, 9]), "sim.bonds correlation dimensions", failures)
      call assert_true(maxval(abs([(bonds_cor(i, i), i = 1, 9)] - 1.0_dp)) < 1.0e-12_dp, &
         "sim.bonds unit diagonal", failures)
      do i = 1, 4
         call assert_close(bonds_loadings(5 + i, i), 0.8_dp, 1.0e-12_dp, "sim.bonds criterion loading", failures)
      end do

      call simulate_vss_data(20, 8, 2, 0.5_dp, vss_items, dichotomous=.true., cut=0.0_dp, info=info)
      call assert_true(info == 0 .and. all(shape(vss_items) == [20, 8]), "sim.VSS dimensions", failures)
      call assert_true(all(abs(vss_items - anint(vss_items)) < 1.0e-12_dp) .and. &
         minval(vss_items) >= 0.0_dp .and. maxval(vss_items) <= 1.0_dp, "sim.VSS dichotomous output", failures)

      call simulate_npl_items(3, [0.0_dp, 1.0_dp], [1.0_dp, 1.0_dp], npl_items, npl_theta, &
         gamma=[1.0_dp], zeta=[1.0_dp], theta_input=[-1.0_dp, 0.0_dp, 1.0_dp], trials=3, info=info)
      call assert_true(info == 0 .and. all(abs(npl_items - 3.0_dp) < 1.0e-12_dp), &
         "sim.npl upper-asymptote responses", failures)
      call assert_true(maxval(abs(npl_theta - [-1.0_dp, 0.0_dp, 1.0_dp])) < 1.0e-12_dp, &
         "sim.npl supplied theta", failures)

      call simulate_npn_items(3, [0.0_dp, 1.0_dp], [1.0_dp, 1.0_dp], npn_items, npn_theta, &
         gamma=[0.0_dp], zeta=[0.0_dp], theta_input=[-1.0_dp, 0.0_dp, 1.0_dp], trials=2, info=info)
      call assert_true(info == 0 .and. maxval(abs(npn_items)) < 1.0e-12_dp, &
         "sim.npn lower-asymptote responses", failures)
      call assert_true(maxval(abs(npn_theta - [-1.0_dp, 0.0_dp, 1.0_dp])) < 1.0e-12_dp, &
         "sim.npn supplied theta", failures)

      latent_cor2 = reshape([1.0_dp, 0.3_dp, 0.3_dp, 1.0_dp], [2, 2])
      marginal_thresholds(1, :) = [-0.5_dp, 0.5_dp]
      marginal_thresholds(2, :) = [-0.5_dp, 0.5_dp]
      call simulate_polychoric_marginals(latent_cor2, marginal_thresholds, 25, poly_marginal, info)
      call assert_true(info == 0 .and. all(shape(poly_marginal) == [25, 2]), "sim.poly.mat dimensions", failures)
      call assert_true(minval(poly_marginal) >= 0.0_dp .and. maxval(poly_marginal) <= 2.0_dp .and. &
         all(abs(poly_marginal - anint(poly_marginal)) < 1.0e-12_dp), "sim.poly.mat categories", failures)

      call block_random_design(8, design, info, conditions=[2, 2])
      call assert_true(info == 0 .and. all(shape(design) == [8, 3]), "block.random dimensions", failures)
      call assert_true(all(design(1:4, 1) == 1) .and. all(design(5:8, 1) == 2), "block.random block labels", failures)
      do i = 0, 1
         duplicates = 0
         do j = 1, 4
            duplicates = duplicates + count(design(i * 4 + 1:i * 4 + 4, 2) == design(i * 4 + j, 2) .and. &
               design(i * 4 + 1:i * 4 + 4, 3) == design(i * 4 + j, 3))
         end do
         call assert_true(duplicates == 4, "block.random unique condition combinations", failures)
      end do

      data_summary(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
      data_summary(:, 2) = [2.0_dp, ieee_value(0.0_dp, ieee_quiet_nan), 4.0_dp, 6.0_dp]
      data_summary(:, 3) = [1.0_dp, 2.0_dp, ieee_value(0.0_dp, ieee_quiet_nan), 4.0_dp]
      call describe_fast_numeric(data_summary, valid_counts, complete_cases)
      call assert_true(all(valid_counts == [4, 3, 3]), "describeFast valid counts", failures)
      call assert_true(complete_cases == 2, "describeFast complete cases", failures)
      call describe_data_numeric(data_summary, valid_counts, complete_cases, head_values, tail_values, &
         head_count=2, tail_count=2)
      call assert_true(all(shape(head_values) == [3, 2]), "describeData head dimensions", failures)
      call assert_close(head_values(1, 2), 2.0_dp, 1.0e-12_dp, "describeData transposed head", failures)
      call assert_close(tail_values(1, 2), 4.0_dp, 1.0e-12_dp, "describeData transposed tail", failures)

      xdata(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
      xdata(:, 2) = [4.0_dp, 3.0_dp, 2.0_dp, 1.0_dp]
      call big_correlation_matrix(xdata, big_r, block_size=1, method="pearson", info=info)
      call assert_true(info == 0 .and. all(shape(big_r) == [2, 2]), "bigCor dimensions", failures)
      call assert_close(big_r(1, 2), -1.0_dp, 1.0e-12_dp, "bigCor Pearson block stitching", failures)
      call big_correlation_matrix(xdata, big_r, block_size=1, method="spearman", info=info)
      call assert_close(big_r(1, 2), -1.0_dp, 1.0e-12_dp, "bigCor Spearman mode", failures)

      do i = 1, 10
         mixed_data(i, 1) = real(i, dp)
         mixed_data(i, 2) = real(1 + mod(i - 1, 3), dp)
         mixed_data(i, 3) = real(merge(0, 1, i <= 5), dp)
      end do
      mixed_types = [1, 2, 3]
      call mixed_correlation_matrix(mixed_data, mixed_r, variable_types=mixed_types, detected_types=detected_types, info=info)
      call assert_true(info == 0 .and. all(detected_types == mixed_types), "mixedCor explicit types", failures)
      call assert_true(all(shape(mixed_r) == [3, 3]), "mixedCor dimensions", failures)
      call assert_close(mixed_r(1, 1), 1.0_dp, 1.0e-12_dp, "mixedCor unit diagonal", failures)
      call assert_true(mixed_r(1, 3) > 0.9_dp .and. mixed_r(2, 3) > -1.0_dp, &
         "mixedCor mixed-family estimates", failures)
      call mixed_correlation_matrix(mixed_data, mixed_r, ncat=3, detected_types=detected_types, info=info)
      call assert_true(all(detected_types == mixed_types), "mixedCor automatic types", failures)

      call polydi_correlations(mixed_data(:, 2:2), mixed_data(:, 3:3), polydi_r, info=info)
      call assert_true(info == 0 .and. all(shape(polydi_r) == [1, 1]), "polydi dimensions", failures)
      call assert_close(polydi_r(1, 1), mixed_r(2, 3), 2.0e-10_dp, "polydi matches mixed categorical block", failures)

      group_data(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
      group_labels = [1, 1, 1, 2, 2, 2]
      call stats_by_numeric(group_data, group_labels, group_values, group_means, group_sds, group_counts, &
         grand_means, ms_between, ms_within, f_values, icc1, icc2, info)
      call assert_true(info == 0 .and. all(group_values == [1, 2]), "statsBy group labels", failures)
      call assert_close(group_means(1, 1), 2.0_dp, 1.0e-12_dp, "statsBy first mean", failures)
      call assert_close(group_means(2, 1), 5.0_dp, 1.0e-12_dp, "statsBy second mean", failures)
      call assert_close(grand_means(1), 3.5_dp, 1.0e-12_dp, "statsBy grand mean", failures)
      call assert_close(ms_between(1), 13.5_dp, 1.0e-12_dp, "statsBy between mean square", failures)
      call assert_close(ms_within(1), 1.0_dp, 1.0e-12_dp, "statsBy within mean square", failures)
      call assert_close(f_values(1), 13.5_dp, 1.0e-12_dp, "statsBy F ratio", failures)
      call assert_close(icc1(1), 12.5_dp / 15.5_dp, 1.0e-12_dp, "statsBy ICC1", failures)
      call assert_close(icc2(1), 12.5_dp / 13.5_dp, 1.0e-12_dp, "statsBy ICC2", failures)
   end subroutine test_additional_translations


   subroutine test_simulation_algorithms(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      real(dp) :: f(2, 1), gload(2), fload(4, 2), explicit_major(4, 2), rw(2, 2), rb(2, 2), eta(2)
      real(dp), allocatable :: model(:, :), reliability(:), sl(:, :), generated(:, :)
      real(dp), allocatable :: sphere(:, :), sphere_loadings(:, :), latent(:, :), within(:, :), between(:, :), obs(:, :)
      real(dp), allocatable :: iv1(:), iv2(:), iv3(:), y(:), probability(:, :), probability2(:, :)
      integer, allocatable :: group(:), poly_items(:, :), poly_items2(:, :)
      integer :: info

      f(:, 1) = [0.8_dp, 0.6_dp]
      call simulate_structural_model(f, model, reliability, info=info)
      call assert_true(info == 0, "sim structural model status", failures)
      call assert_close(model(1, 2), 0.48_dp, 1.0e-12_dp, "sim structural covariance", failures)
      call assert_close(reliability(1), 0.64_dp, 1.0e-12_dp, "sim structural reliability", failures)

      gload = [0.8_dp, 0.6_dp]
      fload = reshape([0.7_dp, 0.6_dp, 0.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.7_dp, 0.6_dp], [4, 2])
      call simulate_hierarchical_data(gload, fload, model, sl, info=info)
      call assert_true(info == 0 .and. all(shape(sl) == [4, 3]), "sim.hierarchical dimensions", failures)
      call assert_close(model(1, 3), 0.2352_dp, 1.0e-12_dp, "sim.hierarchical cross-factor correlation", failures)
      call assert_close(sl(1, 1), 0.56_dp, 1.0e-12_dp, "sim.hierarchical general loading", failures)
      call assert_close(sl(1, 2), 0.42_dp, 1.0e-12_dp, "sim.hierarchical residual group loading", failures)

      explicit_major = reshape([0.8_dp, 0.6_dp, 0.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.7_dp, 0.5_dp], [4, 2])
      call simulate_minor_factors(4, 2, generated, model, big_loadings=explicit_major, n_small=0, info=info)
      call assert_true(info == 0 .and. all(shape(generated) == [4, 2]), "sim.minor explicit major loadings", failures)
      call assert_close(model(1, 2), 0.48_dp, 1.0e-12_dp, "sim.minor within-factor correlation", failures)
      call assert_close(model(1, 3), 0.0_dp, 1.0e-12_dp, "sim.minor between-factor correlation", failures)

      call simulate_spherical_items(sphere, nx=3, ny=5, n_subjects=12, loadings=sphere_loadings, latent=latent, info=info)
      call assert_true(info == 0 .and. all(shape(sphere) == [12, 12]), "sim.spherical response dimensions", failures)
      call assert_true(all(shape(sphere_loadings) == [12, 4]), "sim.spherical loading dimensions", failures)
      call assert_true(all(shape(latent) == [12, 4]), "sim.spherical latent dimensions", failures)

      rw = reshape([1.0_dp, 0.25_dp, 0.25_dp, 1.0_dp], [2, 2])
      rb = reshape([1.0_dp, 0.50_dp, 0.50_dp, 1.0_dp], [2, 2])
      eta = [0.8_dp, 0.6_dp]
      call simulate_multilevel_data(rw, rb, eta, 2, 8, within, between, group, obs, info)
      call assert_true(info == 0 .and. all(shape(obs) == [8, 2]), "sim.multilevel response dimensions", failures)
      call assert_true(all(group(1:4) == 1) .and. all(group(5:8) == 2), "sim.multilevel group replication", failures)
      call assert_close(obs(1, 1), within(1, 1) * 0.64_dp + between(1, 1) * 0.36_dp, 1.0e-12_dp, &
         "sim.multilevel level weighting", failures)

      call simulate_anova_data(1, 2, 2, 1, y, iv1, iv2, iv3, standardize=.false., info=info)
      call assert_true(info == 0 .and. size(y) == 4, "sim.anova factorial size", failures)
      call assert_true(maxval(abs(iv1 - [-1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp])) < 1.0e-12_dp, &
         "sim.anova first contrast pattern", failures)
      call assert_true(maxval(abs(iv2 - [-1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp])) < 1.0e-12_dp, &
         "sim.anova second contrast pattern", failures)
      call assert_true(maxval(abs(iv3)) < 1.0e-12_dp, "sim.anova one-level third factor", failures)

      call simulate_poly_ideal_logistic(3, 2, poly_items, probability, theta=[0.0_dp, 1.0_dp, -1.0_dp], &
         discrimination=[1.0_dp, 1.0_dp], difficulty=[0.0_dp, 1.0_dp], categories=5, info=info)
      call assert_true(info == 0 .and. all(shape(poly_items) == [3, 2]), "sim.poly.ideal.npl dimensions", failures)
      call assert_close(probability(1, 1), 0.5_dp, 1.0e-12_dp, "sim.poly.ideal.npl peak probability", failures)
      call assert_true(all(poly_items >= 0) .and. all(poly_items <= 4), "sim.poly.ideal.npl category bounds", failures)

      call simulate_poly_ideal_normal(3, 2, poly_items2, probability2, theta=[0.0_dp, 1.0_dp, -1.0_dp], &
         discrimination=[1.0_dp, 1.0_dp], difficulty=[0.0_dp, 1.0_dp], categories=5, info=info)
      call assert_true(info == 0 .and. all(shape(poly_items2) == [3, 2]), "sim.poly.ideal.npn dimensions", failures)
      call assert_close(probability2(1, 1), 0.5_dp, 1.0e-12_dp, "sim.poly.ideal.npn center probability", failures)

      call simulate_poly_ideal(3, 2, poly_items2, probability2, normal_model=.false., &
         theta=[0.0_dp, 1.0_dp, -1.0_dp], discrimination=[1.0_dp, 1.0_dp], difficulty=[0.0_dp, 1.0_dp], info=info)
      call assert_close(probability2(1, 1), 0.5_dp, 1.0e-12_dp, "sim.poly.ideal logistic dispatch", failures)
   end subroutine test_simulation_algorithms


   subroutine test_omega_extensions(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      real(dp) :: r(4, 4), expected_omega, omega_h, omega_total, omega_limit, alpha, lambda6, ecv
      real(dp) :: glb_value, alpha_pc, r_smc, ten_berge(4), scale_r, raw_alpha(2), std_alpha(2), g6(2)
      real(dp) :: t1(4, 3), t2(4, 3)
      real(dp), allocatable :: sl(:, :), general(:), groups(:, :), phi(:, :), higher(:), h2(:), u2(:), p2(:)
      real(dp), allocatable :: complexity(:), group_omega(:, :), direct(:, :), target(:, :), communalities(:)
      real(dp), allocatable :: item_r(:), profile_r(:)
      integer :: info, nf

      r = 0.36_dp
      r(1, 1) = 1.0_dp
      r(2, 2) = 1.0_dp
      r(3, 3) = 1.0_dp
      r(4, 4) = 1.0_dp
      expected_omega = (4.0_dp * 0.6_dp) ** 2 / sum(r)

      call schmid_leiman_decomposition(r, 1, sl, general, groups, phi, higher, h2, u2, p2, complexity, info)
      call assert_true(info == 0 .and. all(shape(sl) == [4, 2]), "schmid one-factor dimensions", failures)
      call assert_true(maxval(abs(abs(general) - 0.6_dp)) < 2.0e-4_dp, "schmid recovers equicorrelation loadings", failures)
      call assert_true(maxval(abs(groups)) < 2.0e-8_dp, "schmid one-factor residual group is zero", failures)

      call omega_from_correlation(r, 1, omega_h, omega_total, omega_limit, alpha, lambda6, ecv, sl, &
         group_omega, info)
      call assert_true(info == 0, "omega one-factor status", failures)
      call assert_close(omega_h, expected_omega, 5.0e-4_dp, "omega hierarchical equicorrelation", failures)
      call assert_close(omega_total, expected_omega, 5.0e-4_dp, "omega total equicorrelation", failures)
      call assert_close(alpha, expected_omega, 5.0e-4_dp, "omega alpha equicorrelation", failures)
      call assert_close(ecv, 1.0_dp, 1.0e-10_dp, "omega one-factor ECV", failures)

      call glb_factor_reliability(r, glb_value, communalities, nf, info)
      call assert_true(info == 0 .and. nf == 1, "glb.fa factor count", failures)
      call assert_close(glb_value, expected_omega, 5.0e-4_dp, "glb.fa equicorrelation", failures)
      call glb_reliability_summary(r, glb_value, alpha_pc, r_smc, ten_berge, communalities, nf, info)
      call assert_true(info == 0, "glb summary status", failures)
      call assert_close(alpha_pc, 1.0_dp - 1.0_dp / 2.08_dp, 1.0e-8_dp, "glb alpha-PC", failures)
      call assert_true(all(ten_berge > 0.0_dp), "glb Ten Berge bounds positive", failures)
      call assert_true(r_smc > 0.0_dp .and. r_smc < 1.0_dp, "glb SMC reliability bounded", failures)

      call direct_schmid_leiman(r, 1, direct, target, info, cut=0.1_dp)
      call assert_true(info == 0 .and. all(shape(direct) == [4, 2]), "directSl dimensions", failures)
      call omega_direct_from_correlation(r, 1, omega_h, omega_total, omega_limit, alpha, lambda6, ecv, &
         direct, group_omega, info, cut=0.1_dp)
      call assert_true(info == 0 .and. omega_h > 0.0_dp .and. omega_total > 0.0_dp, &
         "omegaDirect positive reliability", failures)

      t1 = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
         2.0_dp, 1.0_dp, 4.0_dp, 3.0_dp, &
         3.0_dp, 4.0_dp, 1.0_dp, 2.0_dp], [4, 3])
      t2 = 2.0_dp * t1 + 1.0_dp
      call test_retest_summary(t1, t2, scale_r, raw_alpha, std_alpha, g6, item_r, profile_r, info)
      call assert_true(info == 0, "testRetest numeric status", failures)
      call assert_close(scale_r, 1.0_dp, 1.0e-12_dp, "testRetest scale correlation", failures)
      call assert_true(maxval(abs(item_r - 1.0_dp)) < 1.0e-12_dp, "testRetest item correlations", failures)
      call assert_true(maxval(abs(profile_r - 1.0_dp)) < 1.0e-12_dp, "testRetest profile correlations", failures)
      call assert_close(raw_alpha(1), raw_alpha(2), 1.0e-12_dp, "testRetest raw alpha scale invariance", failures)
      call assert_close(std_alpha(1), std_alpha(2), 1.0e-12_dp, "testRetest standardized alpha invariance", failures)
   end subroutine test_omega_extensions


   subroutine test_parallel_extensions(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      real(dp) :: r(6, 6), data(6, 3), objective
      real(dp), allocatable :: observed_pc(:), pc_mean(:), pc_quantile(:), observed_fa(:), fa_mean(:), fa_quantile(:)
      real(dp), allocatable :: loadings(:, :), uniqueness(:), subject_means(:), within_r(:)
      integer, allocatable :: seed(:)
      integer :: i, j, info, ncomp, nfact, nseed

      r = 0.1_dp
      do i = 1, 6
         r(i, i) = 1.0_dp
      end do
      do j = 1, 3
         do i = 1, 3
            if (i /= j) r(i, j) = 0.6_dp
         end do
      end do
      do j = 4, 6
         do i = 4, 6
            if (i /= j) r(i, j) = 0.6_dp
         end do
      end do
      call random_seed(size=nseed)
      allocate(seed(nseed))
      seed = 137
      call random_seed(put=seed)
      call parallel_analysis_normal(r, 400, 20, 0.95_dp, observed_pc, pc_mean, pc_quantile, observed_fa, &
         fa_mean, fa_quantile, ncomp, nfact, info)
      call assert_true(info == 0 .and. size(observed_pc) == 6, "fa.parallel normal-reference dimensions", failures)
      call assert_close(sum(observed_pc), 6.0_dp, 1.0e-10_dp, "fa.parallel observed eigenvalue trace", failures)
      call assert_close(sum(pc_mean), 6.0_dp, 1.0e-10_dp, "fa.parallel simulated mean eigenvalue trace", failures)
      call assert_true(ncomp >= 2 .and. ncomp <= 6, "fa.parallel component count", failures)
      call assert_true(nfact >= 2 .and. nfact <= 6, "fa.parallel factor count", failures)

      data = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, &
         2.0_dp, 1.0_dp, 4.0_dp, 3.0_dp, 6.0_dp, 5.0_dp, &
         4.0_dp, 3.0_dp, 2.0_dp, 1.0_dp, 5.0_dp, 6.0_dp], [6, 3])
      call random_effect_factor(data, 1, loadings, uniqueness, subject_means, within_r, objective, info, &
         add_fix_noise=.false.)
      call assert_true(info == 0 .and. all(shape(loadings) == [3, 1]), "fa.random loading dimensions", failures)
      call assert_close(subject_means(1), 7.0_dp / 3.0_dp, 1.0e-12_dp, "fa.random first subject mean", failures)
      call assert_true(all(abs(within_r) <= 1.0_dp), "fa.random within correlations bounded", failures)
      call assert_true(all(uniqueness >= 0.0_dp), "fa.random nonnegative uniqueness", failures)
   end subroutine test_parallel_extensions


   subroutine test_selection_extensions(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      real(dp) :: r(5, 5), x(8, 2)
      real(dp), allocatable :: values(:, :), keys(:, :), d(:, :), hedges(:, :), r_equiv(:, :), t_value(:, :)
      real(dp), allocatable :: fits(:, :), benchmark(:, :)
      integer, allocatable :: selected(:, :), counts(:), strata_values(:), n_first(:, :), n_second(:, :)
      integer, allocatable :: parallel_results(:, :), seed(:)
      integer :: criteria(1), group(8), strata(8), info, nseed

      r = reshape([ &
         1.0_dp, 0.8_dp, -0.7_dp, 0.2_dp, 0.1_dp, &
         0.8_dp, 1.0_dp, -0.6_dp, 0.1_dp, 0.2_dp, &
         -0.7_dp, -0.6_dp, 1.0_dp, 0.3_dp, 0.1_dp, &
         0.2_dp, 0.1_dp, 0.3_dp, 1.0_dp, 0.9_dp, &
         0.1_dp, 0.2_dp, 0.1_dp, 0.9_dp, 1.0_dp], [5, 5])
      criteria = [1]
      call best_items_from_correlation(r, criteria, 0.15_dp, selected, values, counts, n_item=3, info=info)
      call assert_true(info == 0 .and. counts(1) == 3, "bestItems selected count", failures)
      call assert_true(all(selected(1:3, 1) == [2, 3, 4]), "bestItems ranking", failures)
      call assert_close(values(2, 1), -0.7_dp, 1.0e-12_dp, "bestItems preserves sign", failures)

      call parcel_keys_from_correlation(r(1:4, 1:4), 2, keys, info=info)
      call assert_true(info == 0 .and. size(keys, 2) == 2, "parcels pair dimensions", failures)
      call assert_close(sum(abs(keys(:, 1))), 2.0_dp, 1.0e-12_dp, "parcels first pair size", failures)
      call assert_close(sum(abs(keys(:, 2))), 2.0_dp, 1.0e-12_dp, "parcels second pair size", failures)

      x(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 2.0_dp, 3.0_dp, 5.0_dp, 6.0_dp]
      x(:, 2) = [2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 1.0_dp, 2.0_dp, 4.0_dp, 5.0_dp]
      group = [1, 1, 2, 2, 1, 1, 2, 2]
      strata = [10, 10, 10, 10, 20, 20, 20, 20]
      call cohen_d_by_groups(x, group, strata, strata_values, d, hedges, r_equiv, t_value, n_first, n_second, info)
      call assert_true(info == 0 .and. all(strata_values == [10, 20]), "cohen.d.by strata", failures)
      call assert_true(all(n_first == 2) .and. all(n_second == 2), "cohen.d.by counts", failures)
      call assert_true(all(d > 0.0_dp) .and. all(hedges > 0.0_dp), "cohen.d.by positive effects", failures)
      call assert_true(all(r_equiv > 0.0_dp) .and. all(t_value > 0.0_dp), "cohen.d.by derived effects", failures)

      call random_seed(size=nseed)
      allocate(seed(nseed))
      seed = 24681357
      call random_seed(put=seed)
      call simulate_parallel_recovery(1, [4], [1], [20], 2, parallel_results, info)
      call assert_true(info == 0 .and. all(shape(parallel_results) == [1, 7]), "sim.parallel dimensions", failures)
      call assert_true(all(parallel_results(1, [1, 2, 7]) == [1, 4, 20]), "sim.parallel design columns", failures)

      seed = 13579
      call random_seed(put=seed)
      call vss_parallel_random(30, 4, 2, fits, info)
      call assert_true(info == 0 .and. all(shape(fits) == [2, 2]), "VSS.parallel dimensions", failures)
      call assert_true(all(fits >= -1.0e-8_dp) .and. all(fits <= 1.000001_dp), "VSS.parallel bounded fits", failures)

      seed = 97531
      call random_seed(put=seed)
      call simulate_circumplex_benchmark([30], [8], benchmark, info)
      call assert_true(info == 0 .and. all(shape(benchmark) == [1, 14]), "circ.simulation dimensions", failures)
      call assert_close(benchmark(1, 1), 8.0_dp, 1.0e-12_dp, "circ.simulation variable count", failures)
      call assert_close(benchmark(1, 2), 30.0_dp, 1.0e-12_dp, "circ.simulation sample size", failures)
      call assert_true(.not. any(ieee_is_nan(benchmark)), "circ.simulation finite diagnostics", failures)
   end subroutine test_selection_extensions


   subroutine test_iclust_extensions(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      real(dp) :: r(6, 6), loadings(4, 2)
      real(dp), allocatable :: clusters(:, :), similarity(:), alpha(:), beta(:), sorted_loadings(:, :)
      integer, allocatable :: merge_left(:), merge_right(:), merge_sign(:), merge_size(:)
      integer, allocatable :: cluster_order(:), assignments(:), item_order(:), sorted_assignments(:)
      integer :: i, j, n_merges, info

      r = 0.1_dp
      do i = 1, 6
         r(i, i) = 1.0_dp
      end do
      do j = 1, 3
         do i = 1, 3
            if (i /= j) r(i, j) = 0.7_dp
         end do
      end do
      do j = 4, 6
         do i = 4, 6
            if (i /= j) r(i, j) = 0.7_dp
         end do
      end do
      call iclust_hierarchical(r, 2, clusters, merge_left, merge_right, merge_sign, similarity, alpha, beta, &
         merge_size, n_merges, info)
      call assert_true(info == 0 .and. n_merges == 4, "ICLUST merge count", failures)
      call assert_true(all(shape(clusters) == [6, 2]), "ICLUST final dimensions", failures)
      call assert_true(all(abs(sum(abs(clusters), dim=1) - 3.0_dp) < 1.0e-12_dp), &
         "ICLUST recovers block sizes", failures)
      call assert_true(all(beta(1:n_merges) > 0.0_dp), "ICLUST positive merge beta", failures)

      loadings = reshape([0.8_dp, 0.7_dp, 0.2_dp, 0.1_dp, 0.1_dp, 0.2_dp, 0.9_dp, 0.8_dp], [4, 2])
      call iclust_sort_loadings(loadings, 0.3_dp, cluster_order, assignments, item_order, sorted_loadings, info=info)
      allocate(sorted_assignments(size(assignments)))
      sorted_assignments = assignments
      call sort_integer_test(sorted_assignments)
      call assert_true(info == 0 .and. all(sorted_assignments == [1, 1, 2, 2]), &
         "ICLUST.sort assignments", failures)
      call assert_true(all(item_order >= 1) .and. all(item_order <= 4), "ICLUST.sort item permutation", failures)
   end subroutine test_iclust_extensions

   pure subroutine sort_integer_test(values)
      integer, intent(inout) :: values(:) !! Integer vector sorted in ascending order for deterministic test comparisons.
      integer :: i, j, key

      do i = 2, size(values)
         key = values(i)
         j = i - 1
         do while (j >= 1)
            if (values(j) <= key) exit
            values(j + 1) = values(j)
            j = j - 1
         end do
         values(j + 1) = key
      end do
   end subroutine sort_integer_test

end program test_extended
