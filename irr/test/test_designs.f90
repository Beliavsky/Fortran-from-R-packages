! SPDX-License-Identifier: GPL-2.0-or-later
program test_designs
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
   use r_kinds, only: dp
   use irr_common, only: irr_ok, irr_invalid, irr_degenerate, irr_resource
   use irr_designs, only: n_cohen_kappa, n2_cohen_kappa, rel_inter_intra, rel_inter_intra_result
   implicit none
   type(rel_inter_intra_result) :: result, transformed
   real(dp) :: n, raw, nan, small(7, 6), goniometer(29, 6), reference(7), marginal(3)
   integer :: status, i, j

   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   call n_cohen_kappa(0.5_dp, 0.5_dp, 0.85_dp, 0.7_dp, n, status, raw_n = raw)
   call check_sample(n, status, 117.0_dp, 'binary published equal marginal example')
   call require(n - 1.0_dp < raw .and. raw <= n, 'sample size rounding')
   call n_cohen_kappa(0.3_dp, 0.4_dp, 0.55_dp, 0.25_dp, n, status, twosided = .true.)
   call check_sample(n, status, 77.0_dp, 'binary unequal marginals, two-sided')
   call n_cohen_kappa(0.5_dp, 0.5_dp, -0.2_dp, -0.5_dp, n, status, alpha = 0.01_dp, power = 0.9_dp)
   call check_sample(n, status, 119.0_dp, 'binary negative feasible kappas')
   call n_cohen_kappa(0.5_dp, 0.5_dp, 1.0_dp, 0.7_dp, n, status)
   call check_sample(n, status, 16.0_dp, 'binary perfect alternative')
   call n_cohen_kappa(0.1_dp, 0.9_dp, 0.8_dp, 0.0_dp, n, status)
   call require(status == irr_invalid .and. ieee_is_nan(n), 'infeasible binary kappa rejected')
   call n_cohen_kappa(1.0_dp, 1.0_dp, 0.8_dp, 0.0_dp, n, status)
   call require(status == irr_degenerate, 'constant binary marginals')
   call n_cohen_kappa(0.5_dp, 0.5_dp, 0.8_dp, 0.0_dp, n, status, alpha = 0.0_dp)
   call require(status == irr_invalid, 'invalid significance')
   call n_cohen_kappa(nan, 0.5_dp, 0.8_dp, 0.0_dp, n, status)
   call require(status == irr_invalid, 'nonfinite binary marginal')
   call n_cohen_kappa(0.5_dp, 0.5_dp, 0.8_dp, 0.8_dp, n, status)
   call require(status == irr_invalid, 'identical null and alternative')

   marginal = [0.2_dp, 0.25_dp, 0.55_dp]
   call n2_cohen_kappa(marginal, 0.6_dp, 0.4_dp, n, status, raw_n = raw)
   call check_sample(n, status, 101.0_dp, 'three-category published example')
   call require(n - 1.0_dp < raw .and. raw <= n, 'LP sample size rounding')
   call n2_cohen_kappa(marginal, 0.6_dp, 0.4_dp, n, status, power = 0.9_dp)
   call check_sample(n, status, 136.0_dp, 'three-category 90 percent power')
   call n2_cohen_kappa([0.2_dp, 0.05_dp, 0.2_dp, 0.05_dp, 0.2_dp, 0.3_dp], 0.5_dp, 0.1_dp, n, status)
   call check_sample(n, status, 18.0_dp, 'six-category published example')
   call n2_cohen_kappa([0.5_dp, 0.5_dp], 0.85_dp, 0.7_dp, n, status)
   call check_sample(n, status, 117.0_dp, 'LP agrees with binary formula')
   call n2_cohen_kappa([0.25_dp, 0.25_dp, 0.25_dp, 0.25_dp], 0.7_dp, 0.4_dp, n, status, twosided = .true.)
   call check_sample(n, status, 35.0_dp, 'four uniform categories two-sided')
   call n2_cohen_kappa(spread(0.1_dp, 1, 10), 0.6_dp, 0.3_dp, n, status)
   call check_sample(n, status, 20.0_dp, 'maximum supported ten categories')
   call n2_cohen_kappa([0.0_dp, 0.3_dp, 0.7_dp], 0.7_dp, 0.4_dp, n, status)
   call check_sample(n, status, 59.0_dp, 'zero-probability category')
   call n2_cohen_kappa(marginal(3:1:-1), 0.6_dp, 0.4_dp, n, status)
   call check_sample(n, status, 101.0_dp, 'category permutation invariance')
   call n2_cohen_kappa(spread(1.0_dp / 11.0_dp, 1, 11), 0.6_dp, 0.3_dp, n, status)
   call require(status == irr_resource, 'LP category resource bound')
   call n2_cohen_kappa([1.0_dp], 0.6_dp, 0.3_dp, n, status)
   call require(status == irr_invalid, 'one category invalid')
   call n2_cohen_kappa([0.4_dp, 0.5_dp], 0.6_dp, 0.3_dp, n, status)
   call require(status == irr_invalid, 'marginals must sum to one')
   call n2_cohen_kappa([-0.1_dp, 1.1_dp], 0.6_dp, 0.3_dp, n, status)
   call require(status == irr_invalid, 'negative marginal rejected')
   call n2_cohen_kappa([0.5_dp, 0.5_dp], 1.0_dp, 0.3_dp, n, status)
   call require(status == irr_degenerate, 'zero LP variance rejected as upstream')
   call n2_cohen_kappa([0.9_dp, 0.1_dp], -0.9_dp, 0.3_dp, n, status)
   call require(status == irr_invalid, 'infeasible LP kappa rejected')

   goniometer = reshape(real([ &
      -2,16,5,11,7,-7,18,4,0,0,-3,3,7,-6,1,-13,2,4,-10,8,7,-3,-5,5,0,7,-8,1,-3, &
      0,16,6,10,8,-8,19,5,-3,0,-2,-1,9,-7,1,-14,1,4,-9,9,6,-2,-5,5,-1,6,-8,1,-3, &
      1,15,6,10,6,-8,19,5,-2,-2,-2,1,9,-6,0,-14,0,3,-10,8,7,-4,-7,5,-1,6,-8,2,-3, &
      2,12,4,9,5,-9,17,5,-7,1,-4,-1,4,-8,-2,-12,-1,7,-10,2,8,-5,-6,3,-4,4,-10,1,-5, &
      1,14,4,7,6,-10,17,5,-6,2,-3,-2,4,-10,-2,-12,0,6,-11,8,7,-5,-8,4,-3,4,-11,-1,-4, &
      1,13,4,8,6,-9,17,5,-5,1,-3,1,2,-9,-3,-12,0,4,-10,8,7,-5,-7,4,-4,4,-10,0,-5], dp), [29, 6])
   call rel_inter_intra(goniometer, 2, result)
   call require(result%status == irr_ok, 'goniometer repeated ANOVA status')
   call check_vector(result%rohat, [0.945065704196461_dp, 0.984249114673392_dp, &
      0.961251443214779_dp, 0.983979355948605_dp, 0.986438933579531_dp, 0.982068996692559_dp, &
      0.986206153348779_dp, 0.981762591857534_dp], 2.0e-12_dp, 'goniometer reliabilities')
   call check_vector(result%fs, [17.3735406785844_dp, 14.0484917766797_dp, &
      16.3533224587912_dp, 12.3130898513252_dp], 2.0e-12_dp, 'goniometer F statistics')
   call check_vector(result%pvalue, [2.04033456796537e-11_dp, 0.0_dp, 0.0_dp, &
      1.4432899320127e-15_dp], 1.0e-14_dp, 'goniometer F probabilities')
   call check_vector(result%lowvalue, [0.853592283119829_dp, 0.935356881665683_dp, &
      0.973051501584253_dp, 0.976790844808306_dp, 0.969331047117566_dp], 2.0e-11_dp, 'goniometer lower bounds')
   call check_vector(result%sem, [0.92537659311248_dp, 0.857690027870237_dp, &
      0.988438917815801_dp, 1.72817534941168_dp, 1.43915075608538_dp], 2.0e-12_dp, 'goniometer SEMs')

   do j = 1, 6
      do i = 1, 7
         small(i, j) = 0.4_dp * real(i ** 2, dp) + 0.3_dp * real(j, dp) &
            + real(modulo(i * j + j ** 2, 7), dp) * 0.7_dp
      end do
   end do
   call rel_inter_intra(small, 3, result, rho0inter = 0.2_dp, rho0intra = 0.3_dp, conf_level = 0.9_dp)
   call require(result%status == irr_ok, 'three-rater repeated ANOVA status')
   call require(result%subjects == 7 .and. result%raters == 3 .and. result%repetitions == 2, 'design dimensions')
   call check_vector(result%rohat, [0.949445962376724_dp, 0.95695628079327_dp, &
      0.952779499324907_dp, 0.956805152836495_dp, 0.950790707102804_dp, 0.969530470925734_dp, &
      0.950790707102804_dp, 0.950619048666242_dp, 0.969422080688928_dp, 0.950619048666241_dp], &
      2.0e-12_dp, 'three-rater reliabilities')
   call check_vector(result%fs, [64.3161645610626_dp, 23.9674018013756_dp, &
      20.829385266259_dp, 34.3031790358523_dp, 20.829385266259_dp], 2.0e-12_dp, 'three-rater F statistics')
   call check_vector(result%pvalue, [1.99930692001971e-8_dp, 0.00024776067801846_dp, &
      0.000390932720240267_dp, 7.58328053342971e-5_dp, 0.000390932720240267_dp], &
      1.0e-13_dp, 'three-rater probabilities')
   call check_vector(result%lowvalue, [0.903575250350722_dp, 0.909531833028201_dp, &
      0.880545208167668_dp, 0.863774717431673_dp, 0.915008200923793_dp, 0.863774717431672_dp], &
      2.0e-11_dp, 'three-rater lower bounds')
   call check_vector(result%sem, [1.47252617407411_dp, 1.57955689989313_dp, &
      1.23085336250912_dp, 1.57955689989313_dp, 1.59582789381144_dp, 1.53961575292885_dp], &
      2.0e-12_dp, 'three-rater SEMs')
   call pairwise_anova(small, reference)
   call check_vector(result%mean_squares, reference, 2.0e-12_dp, 'independent pairwise ANOVA decomposition')
   call rel_inter_intra(3.0_dp * small + 100.0_dp, 3, transformed, &
      rho0inter = 0.2_dp, rho0intra = 0.3_dp, conf_level = 0.9_dp)
   call require(transformed%status == irr_ok, 'affine transform status')
   call check_vector(transformed%rohat, result%rohat, 2.0e-12_dp, 'affine invariant reliability')
   call check_vector(transformed%lowvalue, result%lowvalue, 2.0e-11_dp, 'affine invariant bounds')
   call check_vector(transformed%fs, result%fs, 2.0e-12_dp, 'affine invariant F tests')
   call check_vector(transformed%sem, 3.0_dp * result%sem, 2.0e-12_dp, 'SEM scales with scores')
   call check_vector(transformed%mean_squares, 9.0_dp * result%mean_squares, 2.0e-12_dp, 'ANOVA square scaling')
   call rel_inter_intra(small, 1, result)
   call require(result%status == irr_invalid, 'two-way model needs multiple raters')
   call rel_inter_intra(small, 4, result)
   call require(result%status == irr_invalid, 'unbalanced design rejected')
   call rel_inter_intra(small(:, 1:3), 3, result)
   call require(result%status == irr_invalid, 'replication required')
   call rel_inter_intra(small, 1001, result)
   call require(result%status == irr_resource, 'rater resource bound')
   call rel_inter_intra(small * 0.0_dp, 3, result)
   call require(result%status == irr_degenerate, 'constant repeated observations')
   small(1, 1) = nan
   call rel_inter_intra(small, 3, result)
   call require(result%status == irr_invalid, 'missing observations explicitly rejected')
   print '(a)', 'irr design and repeated-measures tests passed'

contains
   subroutine require(condition, label)
      logical, intent(in) :: condition !! Assertion condition that must be true for the test to pass.
      character(len = *), intent(in) :: label !! Diagnostic identifying the failed invariant or numerical comparison.

      if (.not. condition) then
         print '(a)', 'FAILED: ' // label
         error stop 1
      end if
   end subroutine require

   subroutine check_sample(actual, sample_status, expected, label)
      real(dp), intent(in) :: actual !! Computed upward-rounded sample size.
      integer, intent(in) :: sample_status !! Status returned by the sample-size procedure.
      real(dp), intent(in) :: expected !! Exact integer-valued result generated by upstream R.
      character(len = *), intent(in) :: label !! Description of the sample-size fixture.

      if (sample_status /= irr_ok .or. actual /= expected) then
         print *, trim(label), ': status, actual, expected = ', sample_status, actual, expected
         error stop 1
      end if
   end subroutine check_sample

   subroutine check_vector(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual(:) !! Computed numerical vector to compare.
      real(dp), intent(in) :: expected(:) !! Independently generated reference vector with the same shape.
      real(dp), intent(in) :: tolerance !! Allowed absolute error times max(1, absolute reference value).
      character(len = *), intent(in) :: label !! Description printed if the vector comparison fails.

      call require(size(actual) == size(expected), label // ' shape')
      call require(all(abs(actual - expected) <= tolerance * max(1.0_dp, abs(expected))), label)
   end subroutine check_vector

   pure subroutine pairwise_anova(x, ms)
      real(dp), intent(in) :: x(7, 6) !! Independent fixture with seven subjects, three raters, and two repeats.
      real(dp), intent(out) :: ms(7) !! MSS, MSR, interaction mean square, pooled residual MSE, and three rater residual MSEs.
      real(dp) :: subject_average(7), rater_average(3), flattened(42), total_ss, subject_ss, rater_ss, residual_ss
      integer :: a, b, r

      subject_average = sum(x, dim = 2) / 6.0_dp
      do r = 1, 3
         rater_average(r) = sum(x(:, 2 * r - 1:2 * r)) / 14.0_dp
         ms(4 + r) = sum((x(:, 2 * r - 1) - x(:, 2 * r)) ** 2) / 14.0_dp
      end do
      ! Pairwise differences give total centered sum of squares without centering.
      flattened = reshape(x, [42])
      total_ss = 0.0_dp
      do a = 1, 42
         do b = a + 1, 42
            total_ss = total_ss + (flattened(a) - flattened(b)) ** 2 / 42.0_dp
         end do
      end do
      subject_ss = 0.0_dp
      do a = 1, 7
         do b = a + 1, 7
            subject_ss = subject_ss + (subject_average(a) - subject_average(b)) ** 2 * 6.0_dp / 7.0_dp
         end do
      end do
      rater_ss = 0.0_dp
      do a = 1, 3
         do b = a + 1, 3
            rater_ss = rater_ss + (rater_average(a) - rater_average(b)) ** 2 * 14.0_dp / 3.0_dp
         end do
      end do
      residual_ss = 7.0_dp * sum(ms(5:))
      ms(1) = subject_ss / 6.0_dp
      ms(2) = rater_ss / 2.0_dp
      ms(3) = (total_ss - subject_ss - rater_ss - residual_ss) / 12.0_dp
      ms(4) = residual_ss / 21.0_dp
   end subroutine pairwise_anova
end program test_designs
