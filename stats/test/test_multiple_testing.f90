program test_multiple_testing
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_quiet_nan, ieee_value
   use r_stats, only: dp, p_adjust
   use stats_test_assertions, only: assert_true, assert_vector_close
   implicit none
   real(dp), parameter :: p(6) = [0.04_dp, 0.001_dp, 0.03_dp, 0.2_dp, 0.03_dp, 0.8_dp]
   real(dp), allocatable :: result(:)
   real(dp) :: missing

   call check(p_adjust(p, 'hommel'), &
      [0.12_dp, 0.006_dp, 0.09_dp, 0.4_dp, 0.09_dp, 0.8_dp], 'Hommel R')
   call check(p_adjust(p, 'hommel', 10), &
      [0.28_dp, 0.01_dp, 0.21_dp, 1.0_dp, 0.21_dp, 1.0_dp], 'Hommel larger n R')
   call check(p_adjust([0.01_dp, 0.04_dp, 0.03_dp, 0.06_dp], 'hommel'), &
      [0.04_dp, 0.06_dp, 0.06_dp, 0.06_dp], 'Hommel distinguishing case R')
   call check(p_adjust([0.04_dp, 0.01_dp], 'hommel'), [0.04_dp, 0.02_dp], 'Hommel pair R')
   call check(p_adjust([0.3_dp], 'hommel'), [0.3_dp], 'Hommel singleton')
   call check(p_adjust([0.0_dp, 1.0_dp], 'hommel'), [0.0_dp, 1.0_dp], 'Hommel endpoints')

   call check(p_adjust(p), [0.15_dp, 0.006_dp, 0.15_dp, 0.4_dp, 0.15_dp, 0.8_dp], 'Holm')
   call check(p_adjust(p, 'hochberg'), [0.12_dp, 0.006_dp, 0.12_dp, 0.4_dp, 0.12_dp, 0.8_dp], 'Hochberg')
   call check(p_adjust(p, 'bonferroni'), [0.24_dp, 0.006_dp, 0.18_dp, 1.0_dp, 0.18_dp, 1.0_dp], 'Bonferroni')
   call check(p_adjust(p, 'BH'), [0.06_dp, 0.006_dp, 0.06_dp, 0.24_dp, 0.06_dp, 0.8_dp], 'BH')
   call check(p_adjust(p, 'BY'), [0.147_dp, 0.0147_dp, 0.147_dp, 0.588_dp, 0.147_dp, 1.0_dp], 'BY')
   call check(p_adjust(p, 'fdr'), p_adjust(p, 'BH'), 'fdr alias')
   call check(p_adjust(p, 'none'), p, 'none')
   call check(p_adjust(p, n=10), [0.28_dp, 0.01_dp, 0.27_dp, 1.0_dp, 0.27_dp, 1.0_dp], 'Holm n')
   call check(p_adjust(p, 'hochberg', 10), [0.28_dp, 0.01_dp, 0.24_dp, 1.0_dp, 0.24_dp, 1.0_dp], 'Hochberg n')
   call check(p_adjust(p, 'bonferroni', 10), [0.4_dp, 0.01_dp, 0.3_dp, 1.0_dp, 0.3_dp, 1.0_dp], 'Bonferroni n')
   call check(p_adjust(p, 'BH', 10), [0.1_dp, 0.01_dp, 0.1_dp, 0.4_dp, 0.1_dp, 1.0_dp], 'BH n')
   call check(p_adjust(p, 'BY', 10), &
      [0.29289682539682538_dp, 0.029289682539682539_dp, 0.29289682539682538_dp, &
       1.0_dp, 0.29289682539682538_dp, 1.0_dp], 'BY n')
   missing = ieee_value(0.0_dp, ieee_quiet_nan)
   result = p_adjust([0.01_dp, missing, 0.04_dp], 'hommel', 5)
   call assert_true(ieee_is_nan(result(2)), 'Hommel NaN position')
   call check(result([1, 3]), [0.05_dp, 0.16_dp], 'Hommel missing larger n R')
   result = p_adjust([missing, missing], 'hommel')
   call assert_true(all(ieee_is_nan(result)), 'Hommel all missing')
   result = p_adjust([real(dp) ::], 'hommel')
   call assert_true(size(result) == 0, 'Hommel empty')
   result = p_adjust([0.01_dp, missing, 0.04_dp], 'BH')
   call assert_true(ieee_is_nan(result(2)), 'NaN position preserved')
   call check(result([1, 3]), [0.02_dp, 0.04_dp], 'default count excludes NaN')
   result = p_adjust([0.01_dp, missing, 0.04_dp], 'BH', 3)
   call check(result([1, 3]), [0.03_dp, 0.06_dp], 'explicit count includes missing hypothesis')
   call check(p_adjust([0.0_dp, 1.0_dp, 0.05_dp]), [0.0_dp, 1.0_dp, 0.1_dp], 'endpoints')
   result = p_adjust([real(dp) ::])
   call assert_true(size(result) == 0, 'empty input')
   result = p_adjust([missing, missing])
   call assert_true(all(ieee_is_nan(result)), 'all missing')
   call check(p_adjust([0.3_dp]), [0.3_dp], 'singleton')
   print *, 'test_multiple_testing: PASS'

contains

   subroutine check(actual, expected, label)
      !! Checks finite adjusted probabilities against deterministic R reference values.
      real(dp), intent(in) :: actual(:) !! Adjusted probabilities returned by Fortran.
      real(dp), intent(in) :: expected(:) !! Reference probabilities from R.
      character(len=*), intent(in) :: label !! Identification of the reference case.

      call assert_true(all(ieee_is_finite(actual)), label//' finite')
      call assert_vector_close(actual, expected, label, 2.0e-14_dp)
   end subroutine check

end program test_multiple_testing
