program test_strucchange
   use r_kinds, only : dp
   use strucchange, only : ave_f_statistic, breakpoint_confint_result
   use strucchange, only : breakpoint_path_result, breakpoint_result
   use strucchange, only : breakpoint_confidence_intervals
   use strucchange, only : cat_l2_bb_critical_value, cat_l2_bb_pvalue, cat_l2_bb_statistic
   use strucchange, only : compute_breakpoint_path, compute_breakpoints
   use strucchange, only : compute_fstats, efp_pvalue, exp_f_statistic, fstats_pvalue
   use strucchange, only : fstats_result, generalized_fluctuation_process
   use strucchange, only : max_mosum_statistic, mre_critical_value, ols_cusum, ols_mosum
   use strucchange, only : recursive_residuals, segmented_fit, sup_f_statistic, sup_lm_statistic
   implicit none

   integer :: failures
   failures = 0
   call test_regression_family(failures)
   call test_pvalues(failures)
   call test_generalized_process(failures)
   if (failures /= 0) then
      error stop "strucchange tests failed"
   end if
   print '(a)', "All strucchange tests passed."
contains
   subroutine test_regression_family(failures)
      integer, intent(inout) :: failures
      real(dp), parameter :: noise(20) = [ &
         0.20_dp, -0.10_dp, 0.05_dp, -0.15_dp, 0.10_dp, -0.05_dp, &
         0.12_dp, -0.08_dp, 0.03_dp, -0.04_dp, 0.05_dp, -0.07_dp, &
         0.09_dp, -0.11_dp, 0.04_dp, -0.02_dp, 0.06_dp, -0.09_dp, &
         0.08_dp, -0.03_dp ]
      real(dp), parameter :: rr_reference(18) = [ &
         0.18371173_dp, -0.02738613_dp, 0.20554805_dp, 0.00345033_dp, &
         0.14394526_dp, -0.06982972_dp, 0.05324472_dp, -0.01752920_dp, &
         1.96797273_dp, 1.43103406_dp, 1.27344586_dp, 0.83504037_dp, &
         0.80670780_dp, 0.59398009_dp, 0.54952404_dp, 0.30596886_dp, &
         0.39777297_dp, 0.21688730_dp ]
      real(dp), allocatable :: coefficients(:, :), fitted(:), mosum(:, :), process(:, :)
      real(dp), allocatable :: residuals(:), rr(:), segment_rss(:)
      real(dp) :: x(20, 2), y(20)
      type(fstats_result) :: fs
      type(breakpoint_result) :: bp
      type(breakpoint_confint_result) :: ci
      type(breakpoint_path_result) :: path
      integer :: i, info

      do i = 1, 20
         x(i, 1) = 1.0_dp
         x(i, 2) = real(i, dp)
         y(i) = 1.0_dp + 0.4_dp * real(i, dp) + noise(i)
         if (i > 10) y(i) = y(i) + 2.0_dp + 0.3_dp * real(i - 10, dp)
      end do

      call recursive_residuals(x, y, rr, info)
      call check_int("recursive_residuals info", info, 0, failures)
      call check_vector("recursive residuals", rr, rr_reference, 6.0e-8_dp, failures)

      call compute_fstats(x, y, fs, 3, 17)
      call check_int("Fstats info", fs%info, 0, failures)
      call check_int("Fstats maximum breakpoint", fs%breakpoint, 10, failures)
      call check_close("supF", sup_f_statistic(fs), 1017.3506539111903_dp, 2.0e-9_dp, failures)
      call check_close("aveF", ave_f_statistic(fs), 82.37027978839889_dp, 2.0e-9_dp, failures)
      call check_close("expF", exp_f_statistic(fs), 505.96727675449296_dp, 2.0e-9_dp, failures)

      call compute_breakpoints(x, y, 4, 1, bp)
      call check_int("breakpoint info", bp%info, 0, failures)
      call check_int("breakpoint location", bp%breakpoints(1), 10, failures)
      call check_close("breakpoint RSS", bp%rss, 0.1540569696969694_dp, 2.0e-12_dp, failures)

      call compute_breakpoint_path(x, y, 4, 3, path)
      call check_int("breakpoint path info", path%info, 0, failures)
      call check_close("two-break RSS", path%rss(2), 0.11156683982683949_dp, 2.0e-12_dp, failures)
      call check_close("three-break RSS", path%rss(3), 0.10525238095238112_dp, 2.0e-12_dp, failures)

      call segmented_fit(x, y, [10], coefficients, fitted, residuals, segment_rss, info)
      call check_int("segmented fit info", info, 0, failures)
      call check_close("segment 1 intercept", coefficients(1, 1), 1.049333333333333_dp, 2.0e-12_dp, failures)
      call check_close("segment 1 slope", coefficients(2, 1), 0.3924848484848485_dp, 2.0e-12_dp, failures)
      call check_close("segment 2 intercept", coefficients(1, 2), 0.011272727272727_dp, 2.0e-12_dp, failures)
      call check_close("segment 2 slope", coefficients(2, 2), 0.6992727272727273_dp, 2.0e-12_dp, failures)
      call check_close("segmented total RSS", sum(segment_rss), 0.1540569696969694_dp, 2.0e-12_dp, failures)

      call breakpoint_confidence_intervals(x, y, [10], 0.95_dp, ci)
      call check_int("breakpoint confidence info", ci%info, 0, failures)
      if (.not. ci%valid(1)) then
         write (*, '(a)') "FAIL breakpoint confidence interval invalid"
         failures = failures + 1
      else
         call check_int("breakpoint confidence lower", ci%intervals(1, 1), 9, failures)
         call check_int("breakpoint confidence center", ci%intervals(1, 2), 10, failures)
         call check_int("breakpoint confidence upper", ci%intervals(1, 3), 11, failures)
      end if

      call ols_cusum(x, y, process, info)
      call check_int("OLS-CUSUM info", info, 0, failures)
      call check_int("OLS-CUSUM size", size(process, 1), 21, failures)
      call check_close("OLS-CUSUM point 10", process(10, 1), -0.299106895_dp, 2.0e-8_dp, failures)
      call check_close("OLS-CUSUM endpoint", process(21, 1), 0.0_dp, 2.0e-12_dp, failures)

      call ols_mosum(x, y, 0.25_dp, mosum, info)
      call check_int("OLS-MOSUM info", info, 0, failures)
      call check_int("OLS-MOSUM size", size(mosum, 1), 16, failures)
      call check_close("OLS-MOSUM first", mosum(1, 1), 0.77629399_dp, 2.0e-8_dp, failures)
      call check_close("OLS-MOSUM sixth", mosum(6, 1), -1.59089501_dp, 2.0e-8_dp, failures)
   end subroutine test_regression_family

   subroutine test_pvalues(failures)
      integer, intent(inout) :: failures

      call check_close("supF p-value", fstats_pvalue(5.0_dp, "supF", 2, 4.0_dp), &
         0.34455294311656076_dp, 3.0e-13_dp, failures)
      call check_close("aveF p-value", fstats_pvalue(5.0_dp, "aveF", 2, 4.0_dp), &
         0.05810835366808122_dp, 3.0e-13_dp, failures)
      call check_close("expF p-value", fstats_pvalue(2.0_dp, "expF", 2, 4.0_dp), &
         0.15759398022278395_dp, 3.0e-13_dp, failures)
      call check_close("Brownian motion p-value", &
         efp_pvalue(1.2_dp, "Brownian motion", "max", 1), &
         0.005895245603366073_dp, 3.0e-13_dp, failures)
      call check_close("Brownian bridge p-value", &
         efp_pvalue(1.0_dp, "Brownian bridge", "max", 1), &
         0.2699996716773545_dp, 3.0e-12_dp, failures)
      call check_close("categorical chi-square p-value", &
         cat_l2_bb_pvalue(16.0_dp, 1, 2), 6.334248366623988e-5_dp, &
         5.0e-11_dp, failures)
      call check_close("categorical chi-square critical value", &
         cat_l2_bb_critical_value(0.05_dp, 1, 2), 3.8414588206941285_dp, &
         5.0e-11_dp, failures)
      call check_close("monitoring RE critical value", mre_critical_value(0.05_dp), &
         2.7954834829151127_dp, 5.0e-9_dp, failures)
   end subroutine test_pvalues

   subroutine test_generalized_process(failures)
      integer, intent(inout) :: failures
      real(dp), parameter :: root_half = 0.70710678118654752440_dp
      real(dp) :: scores(4, 2), category_process(5, 1)
      real(dp), allocatable :: j12(:, :), process(:, :)
      integer :: info

      scores = 0.0_dp
      scores(1, 1) = 1.0_dp
      scores(2, 1) = -1.0_dp
      scores(3, 2) = 2.0_dp
      scores(4, 2) = -2.0_dp
      call generalized_fluctuation_process(scores, process, j12, info)
      call check_int("gefp info", info, 0, failures)
      call check_close("gefp first component", process(2, 1), root_half, 2.0e-12_dp, failures)
      call check_close("gefp second component", process(4, 2), root_half, 2.0e-12_dp, failures)
      call check_close("gefp bridge endpoint", maxval(abs(process(5, :))), 0.0_dp, 2.0e-12_dp, failures)
      call check_close("supLM statistic", sup_lm_statistic(process, 0.25_dp, 0.75_dp), &
         8.0_dp / 3.0_dp, 2.0e-12_dp, failures)
      call check_close("maxMOSUM statistic", max_mosum_statistic(process, 0.5_dp), &
         root_half, 2.0e-12_dp, failures)

      category_process(:, 1) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
      call check_close("categorical L2 statistic", &
         cat_l2_bb_statistic(category_process, [0.5_dp, 0.5_dp]), &
         16.0_dp, 2.0e-12_dp, failures)
   end subroutine test_generalized_process

   subroutine check_close(name, actual, expected, tolerance, failures)
      character(len = *), intent(in) :: name
      real(dp), intent(in) :: actual, expected, tolerance
      integer, intent(inout) :: failures
      real(dp) :: scale

      scale = max(1.0_dp, abs(expected))
      if (abs(actual - expected) > tolerance * scale) then
         write (*, '(a,2(1x,es24.16))') "FAIL "//trim(name)//":", actual, expected
         failures = failures + 1
      end if
   end subroutine check_close

   subroutine check_vector(name, actual, expected, tolerance, failures)
      character(len = *), intent(in) :: name
      real(dp), intent(in) :: actual(:), expected(:), tolerance
      integer, intent(inout) :: failures

      if (size(actual) /= size(expected)) then
         write (*, '(a,2(1x,i0))') "FAIL "//trim(name)//" size:", size(actual), size(expected)
         failures = failures + 1
         return
      end if
      if (maxval(abs(actual - expected)) > tolerance) then
         write (*, '(a,1x,es24.16)') "FAIL "//trim(name)//" max error:", &
            maxval(abs(actual - expected))
         failures = failures + 1
      end if
   end subroutine check_vector

   subroutine check_int(name, actual, expected, failures)
      character(len = *), intent(in) :: name
      integer, intent(in) :: actual, expected
      integer, intent(inout) :: failures

      if (actual /= expected) then
         write (*, '(a,2(1x,i0))') "FAIL "//trim(name)//":", actual, expected
         failures = failures + 1
      end if
   end subroutine check_int
end program test_strucchange
