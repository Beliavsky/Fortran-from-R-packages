program test_edges
   use matrixtests
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_is_nan, ieee_is_finite, &
      ieee_quiet_nan, ieee_positive_inf
   implicit none
   real(dp) :: na, inf, a(2, 6), b(2, 6), t(6), one(1, 6), pair(1, 3), zeros(1, 0)
   type(test_result), allocatable :: r(:), q(:)
   integer :: i
   na = ieee_value(0.0_dp, ieee_quiet_nan)
   inf = ieee_value(0.0_dp, ieee_positive_inf)
   a(1, :) = [1.0_dp, 2.0_dp, 4.0_dp, 5.0_dp, 8.0_dp, 9.0_dp]
   a(2, :) = [2.0_dp, 3.0_dp, 1.0_dp, 6.0_dp, 4.0_dp, 8.0_dp]
   b = a + 0.7_dp
   r = row_t_welch(a, b)
   q = row_t_welch(a * 1.0e100_dp, b * 1.0e100_dp)
   call check(all(q%status == mt_ok), 'large-scale Welch status')
   call check(maxval(abs(q%df - r%df)) < 1.0e-12_dp, 'large-scale Welch degrees')
   call check(maxval(abs(q%pvalue - r%pvalue)) < 1.0e-12_dp, 'large-scale Welch p')
   r = row_cor_pearson(a, b)
   q = row_cor_pearson(a * 1.0e200_dp, b * 1.0e200_dp)
   call check(all(q%status == mt_ok), 'large-scale correlation status')
   call check(maxval(abs(q%cor - r%cor)) < 1.0e-12_dp, 'large-scale correlation value')
   r = row_jarquebera(a)
   q = row_jarquebera(a * 1.0e200_dp)
   call check(all(q%status == mt_ok), 'large-scale Jarque-Bera status')
   call check(maxval(abs(q%statistic - r%statistic)) < 1.0e-12_dp, 'large-scale Jarque-Bera value')
   q = row_f_var(a * 1.0e200_dp, b * 1.0e200_dp)
   call check(all(q%status == mt_nonfinite), 'overflow F variance explicit status')
   q = row_t_welch(a * 1.0e200_dp, b * 1.0e200_dp)
   call check(all(q%status == mt_nonfinite), 'overflow t variance explicit status')
   r = row_t_equalvar(a, b(1:1, :))
   q = row_t_equalvar(a(2:2, :), b(1:1, :))
   call check(abs(r(2)%pvalue - q(1)%pvalue) < 1.0e-14_dp, 'second-sample broadcast')
   r = row_t_paired(a, b(:, 1:5))
   call check(all(r%status == mt_invalid), 'paired shape validation')
   r = row_t_onesample(zeros)
   call check(r(1)%status == mt_insufficient .and. r(1)%obs == 0, 'empty one-sample')
   r = row_t_welch(a, b, alternative = 'invalid')
   call check(all(r%status == mt_invalid), 'alternative validation')
   r = row_t_welch(a, b, conf_level = na)
   call check(all(ieee_is_nan(r%conf_low)) .and. all(r%status == mt_ok), 'disabled confidence intervals')
   r = row_t_welch(a, b, conf_level = 1.0_dp)
   call check(all(.not. ieee_is_finite(r%conf_low)), 'full confidence level endpoint')
   one = 3.0_dp
   r = row_f_var(one, a(1:1, :))
   call check(r(1)%statistic == 0.0_dp .and. r(1)%pvalue == 0.0_dp, 'zero numerator F variance')
   r = row_f_var(a(1:1, :), one)
   call check(r(1)%statistic == inf .and. r(1)%pvalue == 0.0_dp, 'zero denominator F variance')
   a(1, :) = [1.0_dp, 2.0_dp, na, 4.0_dp, 8.0_dp, 9.0_dp]
   b(1, :) = [na, 1.0_dp, 2.0_dp, 3.0_dp, 5.0_dp, na]
   r = row_t_paired(a(1:1, :), b(1:1, :))
   call check(r(1)%obs_x == 5 .and. r(1)%obs_y == 4 .and. r(1)%obs_paired == 3, 'paired missing counts')
   pair(1, :) = [1.0_dp, 1.0_dp, 3.0_dp]
   q = row_t_onesample(pair)
   call check(abs(r(1)%pvalue - q(1)%pvalue) < 1.0e-14_dp, 'paired missing difference test')
   one(1, :) = [1.0_dp, 2.0_dp, 4.0_dp, 5.0_dp, 8.0_dp, 11.0_dp]
   r = row_oneway_equalvar(one, [1, 1, 1, 2, 2, 2])
   q = row_t_equalvar(one(:, 1:3), one(:, 4:6))
   call check(abs(r(1)%statistic - q(1)%statistic ** 2) < 1.0e-12_dp, 'two-group ANOVA equals squared t')
   r = row_oneway_welch(one, [1, 1, 1, 2, 2, 2])
   q = row_t_welch(one(:, 1:3), one(:, 4:6))
   call check(abs(r(1)%pvalue - q(1)%pvalue) < 1.0e-13_dp, 'two-group Welch equals Welch t')
   r = row_oneway_equalvar(one, [1, 1, 0, 2, 2, 0])
   call check(r(1)%obs_tot == 4 .and. r(1)%obs_groups == 2, 'missing group labels')
   r = row_bartlett(one, [1, 1, 2, 3, 3, 3])
   call check(r(1)%obs_tot == 5 .and. r(1)%obs_groups == 2, 'Bartlett drops singleton groups')
   r = row_bartlett(reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, inf], [1, 5]), [1, 1, 2, 2, 3])
   call check(r(1)%status == mt_ok .and. r(1)%pvalue == 1.0_dp, 'Bartlett ignores dropped infinite singleton')
   r = row_oneway_welch(reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, inf], [1, 5]), [1, 1, 2, 2, 3])
   call check(r(1)%status == mt_ok .and. r(1)%statistic == 8.0_dp, 'Welch ignores dropped infinite singleton')
   r = row_wilcoxon_onesample(reshape([1.0_dp, 2.0_dp, 3.0_dp], [1, 3]))
   call check(r(1)%exact .and. abs(r(1)%pvalue - 0.25_dp) < 1.0e-14_dp, 'signed-rank sign enumeration')
   r = row_wilcoxon_twosample(reshape([1.0_dp, 2.0_dp], [1, 2]), reshape([3.0_dp, 4.0_dp], [1, 2]))
   call check(r(1)%exact .and. abs(r(1)%pvalue - 1.0_dp / 3.0_dp) < 1.0e-14_dp, 'rank-sum label enumeration')
   r = row_wilcoxon_onesample(reshape([0.0_dp, 2.0_dp, 3.0_dp], [1, 3]), exact = .true.)
   call check(.not. r(1)%exact .and. r(1)%corrected .and. r(1)%obs == 2, 'zero disables exact signed rank')
   r = row_wilcoxon_onesample(reshape([0.0_dp, 0.0_dp], [1, 2]))
   call check(.not. r(1)%exact .and. r(1)%corrected .and. r(1)%obs == 0, 'all-zero signed-rank method flags')
   r = row_wilcoxon_onesample(zeros)
   call check(r(1)%exact .and. .not. r(1)%corrected, 'empty signed-rank method flags')
   r = row_wilcoxon_twosample(zeros, one)
   call check(r(1)%exact .and. .not. r(1)%corrected, 'empty rank-sum method flags')
   r = row_wilcoxon_twosample(one, one, exact = .true.)
   call check(.not. r(1)%exact .and. r(1)%corrected, 'ties disable exact rank sum')
   t = [(real(4 * (i - 1), dp), i = 1, 6)]
   one(1, :) = 7.0_dp + 3.0_dp * cos(2.0_dp * acos(-1.0_dp) * (t - 4.0_dp) / 24.0_dp)
   r = row_cosinor(one, t)
   call check(abs(r(1)%mesor - 7.0_dp) < 1.0e-12_dp, 'cosinor mesor')
   call check(abs(r(1)%amplitude - 3.0_dp) < 1.0e-12_dp, 'cosinor amplitude')
   call check(abs(r(1)%acrophase - 4.0_dp) < 1.0e-12_dp, 'cosinor acrophase')
   call check(r(1)%pvalue == 0.0_dp, 'perfect cosinor fit')
   q = row_cosinor(one * 1.0e200_dp, t)
   call check(q(1)%status == mt_nonfinite, 'cosinor overflow explicit status')
   r = row_cosinor(one, 0.0_dp * t)
   call check(r(1)%status == mt_degenerate .and. ieee_is_nan(r(1)%pvalue), 'rank-one cosinor design')
   r = row_cosinor(one, t, period = 0.0_dp)
   call check(r(1)%status == mt_invalid, 'invalid cosinor period')
   print *, 'Independent identities, missingness, limits, and scaled-data tests passed.'
contains
   subroutine check(condition, message)
      logical, intent(in) :: condition !! Whether the independent mathematical or validation assertion holds.
      character(*), intent(in) :: message !! Context printed if the assertion fails.
      if (.not. condition) then
         print *, message
         error stop 'edge test failed'
      end if
   end subroutine check
end program test_edges
