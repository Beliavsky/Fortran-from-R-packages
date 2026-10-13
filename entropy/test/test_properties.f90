! SPDX-License-Identifier: GPL-3.0-or-later
program test_properties
   use entropy_module
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_positive_inf
   implicit none
   real(dp) :: p(4), q(4), table(2, 3), joint(2, 2), prior(2, 2), nan, inf, zero(0), x(7)
   real(dp) :: huge_weights(4), tiny_weights(2), expected, large_table(3, 4)
   integer, allocatable :: bins(:), bins2(:, :)
   integer :: status, checks
   type(test_result) :: test1, test2
   checks = 0
   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   inf = ieee_value(0.0_dp, ieee_positive_inf)
   p = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
   q = [4.0_dp, 3.0_dp, 2.0_dp, 1.0_dp]

   call close (entropy_plugin([1.0_dp, 1.0_dp], 'log2'), 1.0_dp, 'one bit')
   call close (entropy_plugin([0.0_dp, 7.0_dp]), 0.0_dp, 'point mass')
   call close (entropy_plugin(p), entropy_plugin(p * 7.0_dp), 'scale invariant entropy')
   call close (kl_plugin(p, p), 0.0_dp, 'self KL')
   call close (chi2_plugin(p, p), 0.0_dp, 'self chi-squared')
   call close (kl_plugin(p, q), kl_plugin(2.0_dp * p, 9.0_dp * q), 'independent KL rescaling')
   call close (entropy_plugin(p, 'log2') * log(2.0_dp), entropy_plugin(p), 'bits to nats')
   call close (entropy_plugin(p, 'log10') * log(10.0_dp), entropy_plugin(p), 'decimal digits to nats')
   call close (entropy_millermadow(p) - entropy_empirical(p), 3.0_dp / 20.0_dp, 'Miller-Madow correction')
   call close (entropy_chaoshen([10.0_dp, 0.0_dp]), 0.0_dp, 'Chao-Shen point mass')
   call close (entropy_chaoshen([1.0_dp, 1.0_dp]), -log(0.25_dp) / 0.875_dp, 'Chao-Shen singleton correction')
   call check(all(abs(freqs_empirical(p) - p / sum(p)) < 1.0e-15_dp), 'normalized empirical counts')
   call check(all(abs(freqs_dirichlet(p, 0.0_dp) - freqs_empirical(p)) < 1.0e-15_dp), 'zero prior')
   call check(all(abs(freqs_dirichlet([0.0_dp, 0.0_dp], 0.5_dp) - 0.5_dp) < 1.0e-15_dp), 'prior-only estimate')
   call check(all(abs(freqs_shrink(p, 0.0_dp) - freqs_empirical(p)) < 1.0e-15_dp), 'zero shrinkage')
   call check(all(abs(freqs_shrink(p, 1.0_dp) - 0.25_dp) < 1.0e-15_dp), 'full shrinkage')
   call close (shrinkage_intensity([1.0_dp, 1.0_dp, 1.0_dp]), 1.0_dp, 'uniform shrinkage')
   call close (shrinkage_intensity([0.0_dp, 1.0_dp]), 1.0_dp, 'single-observation shrinkage')
   call close (shrinkage_intensity([0.0_dp, 0.0_dp]), 1.0_dp, 'empty-sample intensity mirrors R')
   call check(all(ieee_is_nan(freqs_shrink([0.0_dp, 0.0_dp]))), 'empty sample still lacks frequencies')
   call close (entropy_shrink(p, 0.0_dp), entropy_empirical(p), 'entropy zero shrinkage')
   call close (entropy_shrink(p, 1.0_dp), log(4.0_dp), 'maximum entropy under full shrinkage')
   call close (kl_shrink(p, q, 1.0_dp, 1.0_dp), 0.0_dp, 'KL of uniform targets')
   call close (chi2_shrink(p, q, 1.0_dp, 1.0_dp), 0.0_dp, 'chi-square of uniform targets')

   table = reshape([1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp, 3.0_dp, 6.0_dp], [2, 3])
   call close (mi_plugin(table), 0.0_dp, 'independent table MI')
   call close (chi2indep_plugin(table), 0.0_dp, 'independent table chi-square')
   joint = reshape([4.0_dp, 0.0_dp, 0.0_dp, 4.0_dp], [2, 2])
   call close (mi_plugin(joint, 'log2'), 1.0_dp, 'perfect correlation is one bit')
   call close (chi2indep_plugin(joint), 1.0_dp, 'perfect table chi-square')
   expected = entropy_plugin(sum(joint, dim=1)) + entropy_plugin(sum(joint, dim=2)) &
              - entropy_plugin(reshape(joint, [4]))
   call close (mi_plugin(joint), expected, 'mutual information entropy identity')
   call close (mi_plugin(joint), mi_plugin(transpose(joint)), 'symmetric mutual information')
   prior = reshape([0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp], [2, 2])
   call check(all(abs(freqs_dirichlet(joint, prior) - (joint + prior) / sum(joint + prior)) < 1.0e-15_dp), &
              'matrix-specific Dirichlet prior')
   call check(all(ieee_is_nan(freqs_dirichlet(joint, table))), 'matrix prior shape mismatch')
   joint = reshape([0.0_dp, 2.0_dp, 0.0_dp, 3.0_dp], [2, 2])
   call close (mi_plugin(joint), 0.0_dp, 'zero row allowed for KL-based MI')
   call check(ieee_is_nan(chi2indep_plugin(joint)), 'zero row gives upstream Pearson NaN')
   call close (mi_shrink(joint, 1.0_dp), 0.0_dp, 'full table shrinkage independent')
   call close (chi2indep_shrink(joint, 1.0_dp), 0.0_dp, 'full table chi-square independent')

   test1 = gstat(p, q)
   test2 = gstat(p, q, 'log2')
   call close (test1%pval, test2%pval, 'test p-value invariant under unit')
   call close (test1%stat, test2%stat * log(2.0_dp), 'test statistic converted from bits')
   call check(test1%df == 3, 'goodness-of-fit degrees of freedom')
   test1 = gstatindep(table)
   call close (test1%stat, 0.0_dp, 'independence G statistic')
   call close (test1%pval, 1.0_dp, 'independence p-value')
   call check(test1%df == 2, 'table degrees of freedom')
   large_table = reshape([860125.0_dp, 690375.0_dp, 777000.0_dp, 585868.0_dp, &
                          470244.0_dp, 529248.0_dp, 461027.0_dp, 370041.0_dp, &
                          416472.0_dp, 419741.0_dp, 336903.0_dp, 379176.0_dp], [3, 4])
   test1 = gstatindep(large_table)
   call check(abs(test1%stat) < 1.0e-7_dp, 'large exact-independent table has zero G within roundoff')
   call close (test1%pval, 1.0_dp, 'negative MI roundoff amplified by counts retains p one')
   test1 = chi2stat([2.0_dp, 8.0_dp], [5.0_dp, 5.0_dp])
   call close (test1%stat, 3.6_dp, 'Pearson count-statistic identity')
   test1 = gstat([7.0_dp], [1.0_dp])
   call close (test1%pval, 1.0_dp, 'zero degrees of freedom')
   test1 = gstat([1.0_dp, 0.0_dp], [0.0_dp, 1.0_dp])
   call check(.not. ieee_is_finite(test1%stat), 'infinite G statistic')
   call close (test1%pval, 0.0_dp, 'infinite statistic tail')
   test1 = gstat([0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp])
   call check(ieee_is_nan(test1%stat) .and. ieee_is_nan(test1%pval), 'invalid test sample')
   test1 = gstat([huge(1.0_dp) / 2.0_dp, huge(1.0_dp) / 4.0_dp], [2.0_dp, 1.0_dp])
   call close (test1%stat, 0.0_dp, 'large-total G identity avoids premature overflow')
   test1 = gstatindep(reshape([huge(1.0_dp) / 2.0_dp, huge(1.0_dp) / 4.0_dp], [2, 1]))
   call close (test1%stat, 0.0_dp, 'large-total independence G avoids premature overflow')
   test1 = chi2stat(p, q, 'invalid')
   call check(ieee_is_nan(test1%stat) .and. ieee_is_nan(test1%pval), 'invalid test unit')

   call check(ieee_is_nan(entropy_plugin(zero)), 'empty entropy input')
   call check(ieee_is_nan(entropy_plugin([0.0_dp, 0.0_dp])), 'zero total entropy input')
   call check(ieee_is_nan(entropy_plugin([-1.0_dp, 2.0_dp])), 'negative weights')
   call check(ieee_is_nan(entropy_plugin([nan, 1.0_dp])), 'NaN input')
   call check(ieee_is_nan(entropy_plugin([inf, 1.0_dp])), 'infinite input')
   call check(ieee_is_nan(entropy_plugin(p, 'bogus')), 'invalid unit')
   call check(ieee_is_nan(entropy(p, 'NSB')), 'untranslated NSB not silently approximated')
   call check(ieee_is_nan(entropy(p, 'bogus')), 'invalid method')
   call check(all(ieee_is_nan(freqs(p, 'MM'))), 'MM has no frequencies')
   call check(all(ieee_is_nan(freqs(p, 'CS'))), 'CS has no frequencies')
   call check(all(ieee_is_nan(freqs(p, 'NSB'))), 'NSB has no frequencies')
   call check(all(ieee_is_nan(freqs_dirichlet(p, -1.0_dp))), 'negative pseudocount')
   call check(all(ieee_is_nan(freqs_dirichlet(p, nan))), 'NaN pseudocount')
   call check(all(ieee_is_nan(freqs_dirichlet(p, [1.0_dp]))), 'no silent prior recycling')
   call check(all(ieee_is_nan(freqs_shrink(p, -0.1_dp))), 'invalid negative shrinkage')
   call check(all(ieee_is_nan(freqs_shrink(p, 1.1_dp))), 'invalid excess shrinkage')
   call check(all(ieee_is_nan(freqs_shrink(p, nan))), 'NaN shrinkage')
   call check(ieee_is_nan(kl_plugin(p, [1.0_dp])), 'divergence shape mismatch')
   call check(ieee_is_nan(chi2_plugin(p, [1.0_dp])), 'chi-square shape mismatch')
   call check(ieee_is_nan(chi2_plugin([0.0_dp, 1.0_dp], [0.0_dp, 1.0_dp])), 'common zero Pearson NaN')
   call close (kl_plugin([0.0_dp, 1.0_dp], [0.0_dp, 1.0_dp]), 0.0_dp, 'common zero KL contributes zero')
   call check(.not. ieee_is_finite(kl_plugin([1.0_dp, 0.0_dp], [0.0_dp, 1.0_dp])), 'positive-over-zero KL')
   call check(.not. ieee_is_finite(chi2_plugin([1.0_dp, 0.0_dp], [0.0_dp, 1.0_dp])), 'positive-over-zero chi-square')
   call check(ieee_is_nan(entropy_millermadow([0.0_dp])), 'invalid Miller-Madow total')
   call check(ieee_is_nan(entropy_chaoshen([0.0_dp])), 'invalid Chao-Shen total')
   huge_weights = huge(1.0_dp) / 2.0_dp
   call close (entropy_plugin(huge_weights), log(4.0_dp), 'overflow-safe probability normalization')
   call close (entropy_dirichlet(huge_weights, huge(1.0_dp)), log(4.0_dp), 'overflow-safe prior normalization')
   call check(ieee_is_nan(entropy_millermadow(huge_weights)), 'overflowing sample total rejected')
   tiny_weights = [tiny(1.0_dp), 2.0_dp * tiny(1.0_dp)]
   call close (entropy_plugin(tiny_weights), entropy_plugin([1.0_dp, 2.0_dp]), 'tiny weights normalization')
   call close (kl_plugin([1.0_dp, tiny(1.0_dp)], [tiny(1.0_dp), 1.0_dp]), &
               -log(tiny(1.0_dp)), 'KL computed without ratio overflow')

   x = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp, -1.0_dp, 2.0_dp]
   call discretize(x, 4, bins, [0.0_dp, 1.0_dp], status)
   call check(status == 0 .and. all(bins == [2, 1, 1, 1]), 'right-closed bins with lowest included')
   call discretize([0.0_dp, 0.2_dp, 0.4_dp, 0.6000000000000001_dp, 0.8_dp, 1.0_dp], &
                   5, bins, status=status)
   call check(status == 0 .and. all(bins == [2, 1, 1, 1, 1]), 'R seq rounding at decimal bin boundaries')
   call discretize([0.0_dp, nan, inf, 1.0_dp], 2, bins, [0.0_dp, 1.0_dp], status)
   call check(status == 0 .and. all(bins == [1, 1]), 'explicit range drops nonfinite observations')
   call discretize([0.0_dp, nan, 1.0_dp], 2, bins, status=status)
   call check(status == 1 .and. all(bins == 0), 'default range rejects missing values')
   call discretize([1.0_dp, 1.0_dp], 2, bins, status=status)
   call check(status == 1 .and. all(bins == 0), 'constant sample has invalid breaks')
   call discretize([1.0_dp, nearest(1.0_dp, 1.0_dp)], 2, bins, status=status)
   call check(status == 1, 'indistinguishable adjacent breaks are invalid')
   call discretize(x, 0, bins, status=status)
   call check(status == 1 .and. size(bins) == 0, 'invalid bin count')
   call discretize(x, 2, bins, [1.0_dp, 0.0_dp], status)
   call check(status == 1, 'reversed explicit range')
   call discretize(zero, 2, bins, [0.0_dp, 1.0_dp], status)
   call check(status == 0 .and. all(bins == 0), 'empty sample with valid range')
   call discretize(zero, 2, bins, status=status)
   call check(status == 1, 'empty sample has no default range')
   call discretize([-huge(1.0_dp), 0.0_dp, huge(1.0_dp)], 2, bins, status=status)
   call check(status == 0 .and. all(bins == [2, 1]), 'overflow-safe histogram break interpolation')
   call discretize2d([0.0_dp, 0.5_dp, 1.0_dp], [0.0_dp, 0.5_dp, 1.0_dp], 2, 2, bins2, status=status)
   call check(status == 0 .and. all(bins2 == reshape([2, 0, 0, 1], [2, 2])), 'joint right-closed bins')
   call discretize2d([0.0_dp, 0.5_dp, 1.0_dp], [0.0_dp, nan, 1.0_dp], 2, 2, bins2, &
                     [0.0_dp, 1.0_dp], [0.0_dp, 1.0_dp], status)
   call check(status == 0 .and. sum(bins2) == 2, 'drop pair with missing coordinate')
   call discretize2d([0.0_dp, 1.0_dp], [0.0_dp], 2, 2, bins2, status=status)
   call check(status == 1 .and. all(bins2 == 0), 'joint sample shape mismatch')
   call discretize2d([0.0_dp, 1.0_dp], [0.0_dp, 1.0_dp], -1, 2, bins2, status=status)
   call check(status == 1 .and. size(bins2) == 0, 'joint invalid bins')

   print *, 'Independent identity and edge checks passed:', checks
contains
   subroutine check(condition, label)
      logical, intent(in) :: condition !! Required invariant; a false value terminates the deterministic test.
      character(len=*), intent(in) :: label !! Human-readable name of the invariant or edge-case contract.
      checks = checks + 1
      if (.not. condition) then
         print *, 'FAILED: ', label
         error stop 1
      end if
   end subroutine check

   subroutine close (actual, expected, label)
      real(dp), intent(in) :: actual !! Numerical result produced by the translated API.
      real(dp), intent(in) :: expected !! Independently derived finite reference value.
      character(len=*), intent(in) :: label !! Name of the mathematical identity being checked.
      call check(ieee_is_finite(actual), label)
      call check(abs(actual - expected) < 3.0e-12_dp * max(1.0_dp, abs(expected)), label)
   end subroutine close
end program test_properties
