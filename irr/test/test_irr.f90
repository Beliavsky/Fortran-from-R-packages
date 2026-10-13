! SPDX-License-Identifier: GPL-2.0-or-later
program test_irr
   use irr
   use irr_reference
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
   implicit none
   real(dp) :: x(12, 4), y(30, 4), xx(12, 4, 2), yy(30, 4, 2), alpha_data(4, 30)
   real(dp) :: marginal_fixture(3, 3)
   real(dp) :: missing_x(12, 4), tab(2, 2), alpha_small(3, 3), alpha_extended(3, 4), nan
   real(dp) :: many(1, 1001), many_levels(257, 2), constant(5, 3), perfect(5, 2)
   type(irr_result) :: res, other
   integer :: i, j, failures
   failures = 0
   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   do j = 1, 4
      do i = 1, 12
         x(i, j) = 2.0_dp * sin(real(i, dp) + real(j, dp) / 2.0_dp) + &
            real(i, dp) / 4.0_dp + real(j, dp) / 3.0_dp + cos(real(i * j, dp)) / 3.0_dp
      end do
      do i = 1, 30
         y(i, j) = real(modulo(i * j + i / 3 + j * j, 4) + 1, dp)
      end do
   end do
   do i = 1, 30
      if (modulo(i, 3) /= 0) y(i, 2) = y(i, 1)
   end do
   y(1:4, 1) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
   y(1:4, 2) = y(1:4, 1)
   call check_result('agreement', agree(y), ref_agree(1:1))
   call check_result('agreement tolerance', agree(y, 1.0_dp), ref_agree(2:2))
   call check_result('unweighted kappa', kappa2(y(:, 1:2)), ref_kappa_unweighted)
   call check_result('linear weighted kappa', kappa2(y(:, 1:2), 'equal'), ref_kappa_equal)
   call check_result('squared weighted kappa', kappa2(y(:, 1:2), 'squared'), ref_kappa_squared)
   call check_result('custom weighted kappa', kappa2(y(:, 1:2), costs = [0.0_dp,1.0_dp,2.0_dp,5.0_dp]), ref_kappa_custom)
   call check_result('unsorted custom weights', kappa2(y(:, 1:2), costs = [2.0_dp,0.0_dp,5.0_dp,1.0_dp]), ref_kappa_unsorted)
   call check_result('Fleiss', kappam_fleiss(y), ref_fleiss)
   res = kappam_fleiss(y, detail = .true.)
   call check('category details', real(nint(reshape(res%detail, [12]) * 1000.0_dp), dp) / 1000.0_dp, ref_fleiss_detail)
   call check_result('exact Fleiss', kappam_fleiss(y, exact = .true.), ref_fleiss_exact)
   res = kappam_light(y)
   call check('Light coefficient', [res%value], ref_light)
   call assert('Light negative variance', res%status == irr_degenerate .and. ieee_is_nan(res%p_value))
   call check_result('Maxwell', maxwell(modulo(y(:, 1:2), 2.0_dp)), ref_maxwell)
   call check_result('Bhapkar', bhapkar(y(:, 1:2)), ref_bhapkar)
   call check_result('rater bias', rater_bias_ratings(y(:, 1:2)), ref_rater_bias)
   call check_result('Stuart Maxwell', stuart_maxwell_ratings(y(:, 1:2)), ref_stuart_maxwell)
   call check_result('Kendall', kendall(y), ref_kendall)
   call check_result('Kendall tie correction', kendall(y, .true.), ref_kendall_correct)
   call check_result('mean Pearson Fisher', meancor(x), ref_meancor_true)
   call check_result('mean Pearson raw', meancor(x, .false.), ref_meancor_false)
   call check_result('mean Spearman Fisher', meanrho(y), ref_meanrho_true)
   call check_result('mean Spearman raw', meanrho(y, .false.), ref_meanrho_false)
   call check_result('Finn one-way', finn(x, 7), ref_finn_oneway)
   call check_result('Finn two-way', finn(x, 7, 'twoway'), ref_finn_twoway)
   res = icc(x, 'oneway', 'consistency', 'single', r0 = 0.3_dp, conf_level = 0.8_dp)
   call check_icc('oneway_consistency_single', res, ref_icc_oneway_consistency_single)
   res = icc(x, 'oneway', 'consistency', 'average', r0 = 0.3_dp, conf_level = 0.8_dp)
   call check_icc('oneway_consistency_average', res, ref_icc_oneway_consistency_average)
   res = icc(x, 'oneway', 'agreement', 'single', r0 = 0.3_dp, conf_level = 0.8_dp)
   call check_icc('oneway_agreement_single', res, ref_icc_oneway_agreement_single)
   res = icc(x, 'oneway', 'agreement', 'average', r0 = 0.3_dp, conf_level = 0.8_dp)
   call check_icc('oneway_agreement_average', res, ref_icc_oneway_agreement_average)
   res = icc(x, 'twoway', 'consistency', 'single', r0 = 0.3_dp, conf_level = 0.8_dp)
   call check_icc('twoway_consistency_single', res, ref_icc_twoway_consistency_single)
   res = icc(x, 'twoway', 'consistency', 'average', r0 = 0.3_dp, conf_level = 0.8_dp)
   call check_icc('twoway_consistency_average', res, ref_icc_twoway_consistency_average)
   res = icc(x, 'twoway', 'agreement', 'single', r0 = 0.3_dp, conf_level = 0.8_dp)
   call check_icc('twoway_agreement_single', res, ref_icc_twoway_agreement_single)
   res = icc(x, 'twoway', 'agreement', 'average', r0 = 0.3_dp, conf_level = 0.8_dp)
   call check_icc('twoway_agreement_average', res, ref_icc_twoway_agreement_average)
   call check_result('Robinson', robinson(x), ref_robinson)
   xx(:, :, 1) = x
   xx(:, :, 2) = x + reshape(sin([(real(i, dp), i = 1,48)]), [12, 4])
   call check_result('iota quantitative', iota(xx), ref_iota_quant)
   call check_result('iota standardized', iota(xx, standardize = .true.), ref_iota_stand)
   yy(:, :, 1) = y
   yy(:, :, 2) = modulo(y, 3.0_dp)
   call check_result('iota nominal', iota(yy, 'nominal'), ref_iota_nominal)
   alpha_data = transpose(y)
   call check_result('alpha nominal', kripp_alpha(alpha_data), ref_alpha_nominal)
   call check_result('alpha ordinal', kripp_alpha(alpha_data, 'ordinal'), ref_alpha_ordinal)
   call check_result('alpha interval', kripp_alpha(alpha_data, 'interval'), ref_alpha_interval)
   call check_result('alpha ratio', kripp_alpha(alpha_data, 'ratio'), ref_alpha_ratio)
   alpha_data(1, 1) = nan
   alpha_data(2:4, 2) = nan
   call check_result('alpha missing nominal', kripp_alpha(alpha_data), ref_alpha_missing_nominal)
   call check_result('alpha missing ordinal', kripp_alpha(alpha_data, 'ordinal'), ref_alpha_missing_ordinal)
   call check_result('alpha missing interval', kripp_alpha(alpha_data, 'interval'), ref_alpha_missing_interval)
   call check_result('alpha missing ratio', kripp_alpha(alpha_data, 'ratio'), ref_alpha_missing_ratio)
   missing_x = x
   missing_x(1, 1) = nan
   missing_x(3, 4) = nan
   res = icc(missing_x, 'twoway', 'agreement')
   call check_result('ICC listwise deletion', res, ref_missing_icc)
   call assert('ICC retained count', res%subjects == 10)
   call check_result('correlation listwise deletion', meancor(missing_x), ref_missing_meancor)

   ! Hand-computed binary agreement and McNemar/Bhapkar relationships.
   tab = reshape([20.0_dp, 10.0_dp, 5.0_dp, 15.0_dp], [2, 2])
   call check_result('analytic kappa', kappa2_table(tab), [0.4_dp])
   res = rater_bias(tab)
   call check('analytic bias', [res%value,res%statistic], [1.0_dp / 3.0_dp, 25.0_dp / 15.0_dp])
   call check_result('analytic Stuart Maxwell', stuart_maxwell_mh(tab), [25.0_dp / 15.0_dp])
   call check_result('analytic Bhapkar', bhapkar_table(tab), [25.0_dp / 14.5_dp])
   other = kappa2_table(transpose(tab))
   res = kappa2_table(tab)
   call check('kappa rater symmetry', [other%value,other%statistic], [res%value,res%statistic])
   alpha_small = reshape([1.0_dp,1.0_dp,2.0_dp, 1.0_dp,2.0_dp,2.0_dp, 2.0_dp,2.0_dp,2.0_dp], [3,3])
   res = kripp_alpha(alpha_small, legacy_complete_weights = .false.)
   call check_result('analytic normalized alpha', res, [1.0_dp / 9.0_dp])
   call check('analytic coincidence matrix', reshape(res%coincidence,[4]), [1.0_dp,2.0_dp,2.0_dp,4.0_dp])
   call check_result('legacy complete alpha', kripp_alpha(alpha_small), [1.0_dp / 18.0_dp])
   alpha_extended(:,1:3) = alpha_small
   alpha_extended(:,4) = nan
   call check_result('normalized empty-column invariance', &
      kripp_alpha(alpha_extended, legacy_complete_weights = .false.), [1.0_dp / 9.0_dp])
   marginal_fixture = transpose(reshape([10.0_dp,3.0_dp,2.0_dp, &
      1.0_dp,10.0_dp,4.0_dp,4.0_dp,1.0_dp,10.0_dp], [3,3]))
   call check_result('Stuart Maxwell legacy equal margins', stuart_maxwell_mh(marginal_fixture), [1.8_dp])
   call check_result('Stuart Maxwell corrected equal margins', &
      stuart_maxwell_mh(marginal_fixture, .false.), [5.0_dp / 37.0_dp])
   constant = 3.0_dp
   call check_result('constant agreement', agree(constant), [100.0_dp])
   res = icc(constant)
   call assert('constant ICC degeneracy', res%status == irr_degenerate .and. ieee_is_nan(res%value))
   res = kappa2(constant(:,1:2))
   call assert('constant kappa degeneracy', res%status == irr_degenerate)
   res = meancor(constant)
   call assert('constant correlation degeneracy', res%status == irr_degenerate)
   res = kendall(constant,.true.)
   call assert('all ties corrected W', res%status == irr_degenerate)
   perfect(:,1) = [1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
   perfect(:,2) = perfect(:,1)
   call check_result('perfect agreement kappa', kappa2(perfect), [1.0_dp])
   call check_result('perfect ICC', icc(perfect), [1.0_dp])
   res = meancor(perfect)
   call assert('all perfect correlations dropped', res%status == irr_degenerate .and. res%dropped == 1)
   perfect(:,2) = perfect(:,1) + 1.0_dp
   call check_result('zero exact agreement', agree(perfect), [0.0_dp])
   many = 1.0_dp
   res = agree(many)
   call assert('rater resource limit', res%status == irr_resource)
   many_levels(:,1) = [(real(i,dp),i = 1,257)]
   many_levels(:,2) = many_levels(:,1)
   res = kappa2(many_levels)
   call assert('category resource limit', res%status == irr_resource)
   res = agree(x, -1.0_dp)
   call assert('invalid tolerance', res%status == irr_invalid)
   res = kappa2(y(:,1:2), costs = [0.0_dp,0.0_dp,0.0_dp,0.0_dp])
   call assert('constant costs rejected', res%status == irr_invalid)
   res = icc(x, conf_level = 1.0_dp)
   call assert('invalid confidence', res%status == irr_invalid)
   res = kripp_alpha(reshape([1.0_dp,-1.0_dp],[2,1]), 'ratio')
   call assert('undefined ratio distance', res%status == irr_degenerate)
   if (failures /= 0) error stop 'irr regression failures'
   print *, 'irr agreement and reliability regressions passed'
contains
   subroutine assert(label, condition)
      character(*), intent(in) :: label !! Name identifying a deterministic condition check.
      logical, intent(in) :: condition !! True when the expected contract holds.
      if (.not. condition) then
         print *, 'FAIL: ', label
         failures = failures + 1
      end if
   end subroutine
   subroutine check(label, actual, expected)
      character(*), intent(in) :: label !! Regression case name printed on failure.
      real(dp), intent(in) :: actual(:) !! Computed numerical outputs in the case's documented order.
      real(dp), intent(in) :: expected(:) !! R 0.85 or independent analytic reference values.
      if (any(ieee_is_nan(actual)) .or. any(abs(actual - expected) > 3.0e-10_dp * max(1.0_dp, abs(expected)))) then
         print *, 'FAIL: ', label
         print *, 'actual: ', actual
         print *, 'expected: ', expected
         failures = failures + 1
      end if
   end subroutine
   subroutine check_result(label, actual, expected)
      character(*), intent(in) :: label !! Name of an irr_result numerical regression.
      type(irr_result), intent(in) :: actual !! Computed coefficient and optional inferential results.
      real(dp), intent(in) :: expected(:) !! One value or the value/statistic/p-value triple.
      real(dp) :: output(3)
      call assert(label // ' status', actual%status == irr_ok)
      output = [actual%value, actual%statistic, actual%p_value]
      call check(label, output(1:size(expected)), expected)
   end subroutine
   subroutine check_icc(label, actual, expected)
      character(*), intent(in) :: label !! ICC model/type/unit regression name.
      type(irr_result), intent(in) :: actual !! ICC coefficient, F test, degrees of freedom, and confidence limits.
      real(dp), intent(in) :: expected(:) !! Seven reference values in value/F/p/df1/df2/lower/upper order.
      call assert(label // ' status', actual%status == irr_ok)
      call check(label, [actual%value, actual%statistic, actual%p_value, actual%df1, actual%df2, &
         actual%lbound, actual%ubound], expected)
   end subroutine
end program
