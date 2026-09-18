! SPDX-License-Identifier: MIT
! SPDX-FileComment: Deterministic ANOVA parity tests against R stats.
program test_anova
   use stats_test_assertions, only: assert_close, assert_equal
   use r_stats, only: anova_comparison_t, anova_lm, dp, lm_fit, lm_fit_t
   use r_stats, only: one_way_anova, one_way_anova_t
   implicit none

   integer :: groups(6)
   real(dp) :: full_x(6, 2), reduced_x(6, 1), y(6)
   type(anova_comparison_t) :: comparison
   type(lm_fit_t) :: full, reduced
   type(one_way_anova_t) :: one_way

   reduced_x(:, 1) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
   full_x(:, 1) = reduced_x(:, 1)
   full_x(:, 2) = [1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp]
   y = [-1.8_dp, 4.9_dp, 4.4_dp, 10.2_dp, 10.1_dp, 15.3_dp]
   groups = [1, 1, 2, 2, 3, 3]

   one_way = one_way_anova(y, groups)
   call assert_equal(one_way%status, 0, "one-way ANOVA status")
   call assert_equal(one_way%df_between, 2, "one-way ANOVA numerator df")
   call assert_equal(one_way%df_within, 3, "one-way ANOVA denominator df")
   call assert_close(one_way%ss_between, 124.363333333333_dp, "one-way ANOVA between SS")
   call assert_close(one_way%ss_within, 52.785_dp, "one-way ANOVA within SS")
   call assert_close(one_way%f_statistic, 3.53405323482050_dp, "one-way ANOVA F")
   call assert_close(one_way%p_value, 0.162652290235247_dp, "one-way ANOVA p value")

   reduced = lm_fit(y, reduced_x)
   full = lm_fit(y, full_x)
   comparison = anova_lm(reduced, full)
   call assert_equal(comparison%status, 0, "nested-model ANOVA status")
   call assert_equal(comparison%df_numerator, 1, "nested-model ANOVA numerator df")
   call assert_equal(comparison%df_denominator, 3, "nested-model ANOVA denominator df")
   call assert_close(comparison%sum_of_squares, 13.2859285714286_dp, &
                     "nested-model ANOVA sum of squares")
   call assert_close(comparison%f_statistic, 65.2514909374400_dp, "nested-model ANOVA F")
   call assert_close(comparison%p_value, 0.00396395693351309_dp, "nested-model ANOVA p value")

   print *, "test_anova: PASS"

end program test_anova
