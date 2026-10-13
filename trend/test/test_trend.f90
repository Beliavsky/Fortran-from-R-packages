! SPDX-License-Identifier: GPL-3.0-only
program test_trend
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_nan, ieee_is_finite
   use trend_module
   use trend_r_reference
   implicit none
   type(trend_result) :: a, b
   real(dp) :: matrix(12, 2), null(12, 64), missing, inf
   real(dp), allocatable :: augmented(:)
   integer :: checks, i
   checks = 0
   missing = ieee_value(0.0_dp, ieee_quiet_nan)
   inf = ieee_value(0.0_dp, ieee_positive_inf)
   null = reference_null()
   matrix(:, 1) = rx
   matrix(:, 2) = ry

   call mk_test(rx, a)
   call check_mk('MK two-sided', a, ref_mk_two_sided)
   call mk_test(rx, a, alternative = 'greater')
   call check_mk('MK greater', a, ref_mk_greater)
   call mk_test(rx, a, alternative = 'less')
   call check_mk('MK less', a, ref_mk_less)
   call mk_test(rtied, a, continuity = .false.)
   call check_mk('MK ties no continuity', a, ref_mk_ties)
   call sens_slope(rx, a)
   call close_vector('Sen slope', [a%statistic, a%p_value, a%estimate, a%conf_int], ref_slope)
   call itam(rx, a, conf_level = 0.9_dp)
   call close_vector('ITAM even', [a%statistic, a%p_value, a%estimate, a%standard_error, a%conf_int, &
      reshape(a%sorted_halves, [12])], ref_itam_even)
   call itam(rx * 1.0e200_dp, a, conf_level = 0.9_dp)
   call close_vector('ITAM huge-scale invariance', [a%statistic, a%p_value, a%estimate / 1.0e200_dp, &
      a%standard_error / 1.0e200_dp, a%conf_int / 1.0e200_dp], ref_itam_even(:6))
   call itam(rx(:11), a, conf_level = 0.9_dp)
   call close_vector('ITAM odd midpoint', [a%statistic, a%p_value, a%estimate, a%standard_error, a%conf_int, &
      reshape(a%sorted_halves, [10])], ref_itam_odd)
   call itam([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], a)
   call assert_true(a%p_value == 0.0_dp .and. .not. ieee_is_finite(a%statistic), 'ITAM zero standard error')
   call itam(spread(1.0_dp, 1, 8), a)
   call assert_true(a%status == trend_degenerate, 'ITAM constant halves')
   call itam(rx, a, conf_level = missing)
   call assert_true(a%status == trend_invalid, 'ITAM invalid confidence')
   call sea_sens_slope(rseason, 2, a)
   call close_vector('seasonal Sen', [a%estimate], ref_sea_slope)
   call smk_test(rseason, 2, a)
   call check_seasonal('seasonal MK', a, ref_smk)
   call smk_test(rseason(:23), 2, a, alternative = 'greater', continuity = .false.)
   call check_seasonal('incomplete seasonal MK', a, ref_smk_incomplete)
   call mult_mk_test(matrix, a)
   call close_vector('multisite MK', [a%statistic, a%p_value, a%score, a%variance, &
      reshape(a%covariance, [4])], ref_mult)
   call mult_mk_test(matrix(:8, :), a)
   call close_vector('short multisite MK', [a%statistic, a%p_value, a%score, a%variance, &
      reshape(a%covariance, [4])], ref_mult_short)
   call csmk_test(rseason, 2, a)
   call close_vector('correlated seasonal MK', [a%statistic, a%p_value, a%score, a%variance, &
      reshape(a%covariance, [4])], ref_csmk)
   call partial_mk_test(rx, ry, a)
   call check_partial('partial MK', a, ref_partial_mk_two_sided)
   call partial_mk_test(rx, ry, a, alternative = 'less')
   call check_partial('partial MK less', a, ref_partial_mk_less)
   call partial_mk_test(rx, ry, a, alternative = 'greater')
   call check_partial('partial MK greater', a, ref_partial_mk_greater)
   call partial_cor_trend_test(rtied, ry, a)
   call close_vector('partial Pearson', [a%statistic, a%p_value, a%estimate, &
      reshape(a%covariance, [9])], ref_partial_cor_pearson)
   call partial_cor_trend_test(rtied * 1.0e200_dp, ry * 1.0e200_dp, a)
   call close_vector('partial Pearson huge-scale invariance', [a%statistic, a%p_value, a%estimate, &
      reshape(a%covariance, [9])], ref_partial_cor_pearson)
   call partial_cor_trend_test(rtied, ry, a, method = 'spearman')
   call close_vector('partial Spearman', [a%statistic, a%p_value, a%estimate, &
      reshape(a%covariance, [9])], ref_partial_cor_spearman)
   call bartels_test(rx, a)
   call check_base('Bartels beta', a, ref_bartels)
   call bartels_test(rlong, a)
   call check_base('Bartels normal', a, ref_bartels_long)
   call cs_test(rx, a)
   call check_base('Cox Stuart', a, ref_cs)
   call cs_test(rlong, a)
   call check_base('Cox Stuart large', a, ref_cs_long)
   call wm_test(rx, a)
   call check_base('Wallis Moore', a, ref_wm)
   call wm_test(rlong, a)
   call check_base('Wallis Moore large', a, ref_wm_long)
   call ww_test(rx, a)
   call check_base('Wald Wolfowitz', a, ref_ww)
   call ww_test(rx + 1.0e8_dp, a)
   call check_base('Wald Wolfowitz large-offset invariance', a, ref_ww)
   call ww_test(rx * 1.0e200_dp, a)
   call check_base('Wald Wolfowitz huge-scale invariance', a, ref_ww)
   call ww_test(rlong, a)
   call check_base('Wald Wolfowitz large', a, ref_ww_long)
   call rrod_test(rx(:6), ry(7:), a)
   call check_base('robust rank order', a, ref_rrod_two_sided)
   call rrod_test(rx(:6), ry(7:), a, alternative = 'less')
   call check_base('robust rank order less', a, ref_rrod_less)
   call rrod_test(rx(:6), ry(7:), a, alternative = 'greater')
   call check_base('robust rank order greater', a, ref_rrod_greater)
   call pettitt_test(rtied, a)
   call check_change('Pettitt', a, ref_pettitt)
   call lanzante_test(rx, a)
   call check_change('Lanzante Wilcoxon exact', a, ref_lanzante_wilcox)
   call assert_true(a%exact, 'Lanzante selects exact for small untied samples')
   call lanzante_test(rtied, a)
   call check_change('Lanzante Wilcoxon tied approximation', a, ref_lanzante_ties)
   call assert_true(.not. a%exact, 'Lanzante selects approximation when ties exist')
   call lanzante_test(rx, a, method = 'rrod.test')
   call check_change('Lanzante robust', a, ref_lanzante_rrod)

   ! Identical R-generated normal matrices yield exact Monte Carlo exceedance counts.
   call br_test(rx, a, null_samples = null)
   call check_change('Buishand range R RNG', a, ref_br)
   call bu_test(rx, a, null_samples = null)
   call check_change('Buishand U R RNG', a, ref_bu)
   call snh_test(rx, a, null_samples = null)
   call check_change('SNHT R RNG', a, ref_snh)
   call assert_true(a%simulations == 64, 'supplied null draws counted')
   ! Independent local generator: compare within conservative combined Monte Carlo error.
   call br_test(rx, a, m = 50000, seed = 7781)
   call close_vector('Buishand range independent RNG', [a%p_value], ref_br_large_mc, 0.015_dp)
   call bu_test(rx, a, m = 50000, seed = 7781)
   call close_vector('Buishand U independent RNG', [a%p_value], ref_bu_large_mc, 0.015_dp)
   call snh_test(rx, a, m = 50000, seed = 7781)
   call close_vector('SNHT independent RNG', [a%p_value], ref_snh_large_mc, 0.015_dp)
   call br_test(rx, a, m = 1000, seed = 9901)
   call br_test(rx, b, m = 1000, seed = 9901)
   call assert_true(a%p_value == b%p_value, 'deterministic local RNG')
   call bu_test(rx, a, m = 1000)
   call bu_test(rx, b, m = 1000)
   call assert_true(a%p_value == b%p_value, 'deterministic default seed')

   ! Independent elementary oracles and invariances.
   call mk_test([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], a, continuity = .false.)
   call close_vector('monotone score and variance', [a%score, a%variance, a%tau], &
      [6.0_dp, 26.0_dp / 3.0_dp, 1.0_dp])
   call sens_slope([(2.0_dp * real(i, dp) + 8.0_dp, i = 1, 20)], a)
   call close_vector('exact linear Sen slope', [a%estimate, a%conf_int], [2.0_dp, 2.0_dp, 2.0_dp])
   call mk_test(rx, a)
   call mk_test(-rx, b)
   call close_vector('MK reversal', [a%score, a%p_value], [-b%score, b%p_value])
   call br_test(7.0_dp * rx + 100.0_dp, b, null_samples = null)
   call br_test(rx, a, null_samples = null)
   call close_vector('range affine invariance', [a%statistic, a%p_value], [b%statistic, b%p_value])
   call bu_test(rx + 1.0e8_dp, a, null_samples = null)
   call check_change('Buishand U large-offset invariance', a, ref_bu)
   call snh_test(rx * 1.0e200_dp, a, null_samples = null)
   call check_change('SNHT huge-scale invariance', a, ref_snh)
   call pettitt_test([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], a)
   call assert_true(size(a%change_points) == 4, 'Pettitt returns every maximizing index')
   call assert_true(a%p_value == 1.0_dp, 'constant Pettitt p-value')
   call lanzante_test([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], a)
   call assert_true(a%status == trend_ambiguous, 'Lanzante multiple-split diagnostic')
   call rrod_test([1.0_dp, 2.0_dp], [3.0_dp, 4.0_dp], a)
   call assert_true(a%status == trend_ok .and. a%p_value == 0.0_dp, 'RROD complete separation')
   call assert_true(.not. ieee_is_finite(a%statistic), 'RROD separated z is infinite')

   ! NaN policy differs deliberately between complete-case and omission APIs.
   augmented = [rx(:5), missing, rx(6:)]
   call cs_test(augmented, a)
   call check_base('Cox Stuart missing omission', a, ref_cs)
   call wm_test(augmented, a)
   call check_base('Wallis Moore missing omission', a, ref_wm)
   call ww_test(augmented, a)
   call check_base('Wald Wolfowitz missing omission', a, ref_ww)
   call mk_test(augmented, a)
   call assert_true(a%status == trend_nonfinite, 'MK rejects NaN')
   call rrod_test([rx(:6), missing, inf], [ry(7:), -inf], a)
   call check_base('RROD omits all nonfinite values', a, ref_rrod_two_sided)
   call mk_test([1.0_dp, 1.0_dp, 1.0_dp], a)
   call assert_true(a%status == trend_degenerate .and. ieee_is_nan(a%p_value), 'constant MK diagnosed')
   call mk_test(rx, a, alternative = 'bogus')
   call assert_true(a%status == trend_invalid, 'invalid MK alternative')
   call mk_test([1.0_dp, 2.0_dp], a)
   call assert_true(a%status == trend_insufficient, 'short MK rejected')
   call sens_slope(rx, a, conf_level = 1.0_dp)
   call assert_true(a%status == trend_invalid, 'invalid confidence level')
   call sea_sens_slope(rseason(:23), 2, a)
   call assert_true(a%status == trend_invalid, 'seasonal slope incomplete cycle rejected')
   call smk_test(rseason, 0, a)
   call assert_true(a%status == trend_invalid, 'zero period rejected before division')
   call csmk_test(rseason(:23), 2, a)
   call assert_true(a%status == trend_invalid, 'correlated seasonal incomplete cycle rejected')
   call partial_mk_test(rx, ry(:11), a)
   call assert_true(a%status == trend_invalid, 'partial MK length mismatch')
   call partial_cor_trend_test(rx, ry, a, method = 'invalid')
   call assert_true(a%status == trend_invalid, 'partial correlation invalid method')
   call bartels_test(rx(:9), a)
   call assert_true(a%status == trend_insufficient, 'Bartels needs ten observations')
   call bartels_test(spread(1.0_dp, 1, 12), a)
   call assert_true(a%status == trend_degenerate, 'constant Bartels diagnosed')
   call cs_test([1.0_dp, 1.0_dp, 1.0_dp], a)
   call assert_true(a%status == trend_degenerate, 'Cox Stuart all tied comparisons')
   call ww_test([1.0_dp, 1.0_dp, 1.0_dp], a)
   call assert_true(a%status == trend_degenerate, 'constant Wald Wolfowitz diagnosed')
   call br_test(rx, a, m = 0)
   call assert_true(a%status == trend_invalid, 'zero simulations rejected')
   call bu_test(rx, a, m = 63, null_samples = null)
   call assert_true(a%status == trend_invalid, 'null matrix count mismatch')
   call snh_test(rx, a, null_samples = null(:11, :))
   call assert_true(a%status == trend_invalid, 'null matrix row mismatch')
   call br_test([1.0_dp, 1.0_dp, 1.0_dp], a, m = 10)
   call assert_true(a%status == trend_degenerate, 'constant observed Monte Carlo input')
   null(:, 1) = 1.0_dp
   call bu_test(rx, a, null_samples = null)
   call assert_true(a%status == trend_degenerate, 'constant null realization rejected')
   print '(a,i0,a)', 'trend: ', checks, ' assertions passed'
contains
   subroutine assert_true(condition, name)
      logical, intent(in) :: condition !! Required invariant whose failure terminates the deterministic test program.
      character(len=*), intent(in) :: name !! Human-readable assertion label printed on failure.
      checks = checks + 1
      if (.not. condition) then
         print '(a)', 'FAIL: ' // name
         error stop 1
      end if
   end subroutine assert_true

   subroutine close_vector(name, actual, expected, tolerance)
      character(len=*), intent(in) :: name !! Comparison label identifying the R reference or independent oracle.
      real(dp), intent(in) :: actual(:) !! Computed finite numerical values.
      real(dp), intent(in) :: expected(:) !! Same-shaped reference values from the actual R package or an independent oracle.
      real(dp), intent(in), optional :: tolerance !! Relative-plus-absolute tolerance, default 2e-10; MC comparisons use 0.015.
      real(dp) :: tol
      integer :: k
      tol = 2.0e-10_dp
      if (present(tolerance)) tol = tolerance
      call assert_true(size(actual) == size(expected), name // ' shape')
      do k = 1, size(actual)
         checks = checks + 1
         if (.not. ieee_is_finite(actual(k)) .or. abs(actual(k) - expected(k)) > tol * (1.0_dp + abs(expected(k)))) then
            print '(a,a,i0,2es24.15)', 'FAIL: ', name, k, actual(k), expected(k)
            error stop 1
         end if
      end do
   end subroutine close_vector

   subroutine check_base(name, value, reference)
      character(len=*), intent(in) :: name !! Label identifying the translated test and branch.
      type(trend_result), intent(in) :: value !! Test result expected to report successful finite statistics.
      real(dp), intent(in) :: reference(:) !! R statistic and p-value, in that order.
      call assert_true(value%status == trend_ok, name // ' status')
      call close_vector(name, [value%statistic, value%p_value], reference)
   end subroutine check_base

   subroutine check_mk(name, value, reference)
      character(len=*), intent(in) :: name !! Label for this Mann-Kendall branch.
      type(trend_result), intent(in) :: value !! MK result with score, variance, and tau-b.
      real(dp), intent(in) :: reference(:) !! R z, p-value, score, variance, and tau-b.
      call assert_true(value%status == trend_ok, name // ' status')
      call close_vector(name, [value%statistic, value%p_value, value%score, value%variance, value%tau], reference)
   end subroutine check_mk

   subroutine check_partial(name, value, reference)
      character(len=*), intent(in) :: name !! Label for the partial MK alternative.
      type(trend_result), intent(in) :: value !! Partial MK result with adjusted score and score correlation.
      real(dp), intent(in) :: reference(:) !! R z, p-value, adjusted score, variance, and correlation.
      call assert_true(value%status == trend_ok, name // ' status')
      call close_vector(name, [value%statistic, value%p_value, value%score, value%variance, value%correlation], reference)
   end subroutine check_partial

   subroutine check_seasonal(name, value, reference)
      character(len=*), intent(in) :: name !! Label for this seasonal MK reference case.
      type(trend_result), intent(in) :: value !! Overall and seasonwise numerical test outputs.
      real(dp), intent(in) :: reference(:) !! R overall z,p,S,var followed by seasonwise S,var,z,p,tau vectors.
      call assert_true(value%status == trend_ok, name // ' status')
      call close_vector(name, [value%statistic, value%p_value, value%score, value%variance, value%season_score, &
         value%season_variance, value%season_z, value%season_p, value%season_tau], reference)
   end subroutine check_seasonal

   subroutine check_change(name, value, reference)
      character(len=*), intent(in) :: name !! Label identifying the change-point algorithm and significance method.
      type(trend_result), intent(in) :: value !! Change-point statistic, all maximizing splits, and process values.
      real(dp), intent(in) :: reference(:) !! R statistic, p-value, maximizing splits, and process values in that order.
      call assert_true(value%status == trend_ok, name // ' status')
      call close_vector(name, [value%statistic, value%p_value, real(value%change_points, dp), value%data], reference)
   end subroutine check_change
end program test_trend
