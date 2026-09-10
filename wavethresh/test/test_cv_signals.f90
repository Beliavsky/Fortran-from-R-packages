! SPDX-License-Identifier: GPL-2.0-or-later
program test_cv_signals
   use wavethresh
   implicit none
   real(dp) :: noisy(32)
   real(dp), allocatable :: samples1(:)
   real(dp), allocatable :: samples2(:)
   real(dp), allocatable :: lsw_one(:)
   real(dp), allocatable :: lsw_two(:)
   type(rss_result_t) :: rss
   type(wavelet_cv_result_t) :: cv
   type(wavelet_cv_result_t) :: compiled_cv
   type(wavelet_cv_result_t) :: stationary_cv
   type(wavelet_cv_result_t) :: stationary_level_cv
   type(wd_t) :: spectrum
   type(wd_t) :: checked_one
   type(wd_t) :: checked_two
   type(signal_set_t) :: signals
   real(dp) :: wst_rss
   real(dp) :: wst_level_rss
   integer :: i
   integer :: level

   do i = 1, 32
      noisy(i) = sin(0.21_dp * real(i, dp)) + 0.15_dp * cos(1.17_dp * real(i, dp)) + &
         0.02_dp * real(mod(7 * i, 5) - 2, dp)
   end do
   rss = rsswav(noisy, value=0.2_dp, filter_number=1.0_dp, family="DaubExPhase", threshold_type="soft", ll=1)
   call assert_true(rss%ok, "rsswav status")
   call assert_true(rss%ssq >= 0.0_dp, "rsswav nonnegative")
   cv = wavelet_cv(noisy, filter_number=1.0_dp, family="DaubExPhase", threshold_type="soft", tolerance=0.02_dp, ll=1)
   call assert_true(cv%ok, "WaveletCV status")
   compiled_cv = cwcv(noisy, ll=1, filter_number=1.0_dp, family="DaubExPhase", &
      threshold_type="soft", tolerance=0.02_dp)
   call assert_true(compiled_cv%ok, "CWCV status")
   call assert_true(abs(compiled_cv%cv_threshold - cv%cv_threshold) < 1.0e-12_dp, "CWCV threshold parity")
   call assert_true(maxval(abs(compiled_cv%cv_reconstruction - cv%cv_reconstruction)) < 1.0e-12_dp, &
      "CWCV reconstruction parity")
   call assert_true(cv%cv_threshold >= 0.0_dp, "WaveletCV threshold")
   call assert_true(size(cv%cv_reconstruction) == size(noisy), "WaveletCV reconstruction size")
   wst_rss = get_rss_wst(noisy, 0.2_dp, [1, 2], 1.0_dp, "DaubExPhase", "soft")
   call assert_true(wst_rss >= 0.0_dp .and. wst_rss < huge(1.0_dp), "GetRSSWST objective")
   call assert_true(abs(wvcvlrss(0.2_dp, noisy, [1, 2], 1.0_dp, "DaubExPhase", "soft") - wst_rss) < &
      1.0e-12_dp, "wvcvlrss wrapper")
   stationary_cv = wst_cv(noisy, ll=1, threshold_type="soft", filter_number=1.0_dp, &
      family="DaubExPhase", tolerance=0.02_dp)
   call assert_true(stationary_cv%ok, "wstCV status")
   call assert_true(stationary_cv%cv_threshold >= 0.0_dp, "wstCV threshold")
   call assert_true(size(stationary_cv%cv_reconstruction) == size(noisy), "wstCV reconstruction size")
   wst_level_rss = get_rss_wst_levels(noisy, [0.2_dp, 0.3_dp], [1, 2], 1.0_dp, "DaubExPhase", "soft")
   call assert_true(wst_level_rss >= 0.0_dp .and. wst_level_rss < huge(1.0_dp), "vector GetRSSWST objective")
   call assert_true(abs(wvcvlrss_levels([0.2_dp, 0.3_dp], noisy, [1, 2], 1.0_dp, "DaubExPhase", "soft") - &
      wst_level_rss) < 1.0e-12_dp, "vector wvcvlrss wrapper")
   stationary_level_cv = wst_cvl(noisy, ll=1, threshold_type="soft", filter_number=1.0_dp, &
      family="DaubExPhase", tolerance=0.05_dp)
   call assert_true(stationary_level_cv%ok, "wstCVl status")
   call assert_true(size(stationary_level_cv%cv_thresholds) == 4, "wstCVl threshold vector")
   call assert_true(all(stationary_level_cv%cv_thresholds >= 0.0_dp), "wstCVl nonnegative thresholds")
   call assert_true(size(stationary_level_cv%cv_reconstruction) == size(noisy), "wstCVl reconstruction size")
   call assert_true(abs(full_wavelet_cv(noisy, 1, "soft", 1.0_dp, "DaubExPhase", 0.02_dp) - &
      cv%cv_threshold) < 1.0e-12_dp, "FullWaveletCV threshold wrapper")

   signals = dj_ex(64)
   call assert_true(size(signals%blocks) == 64, "DJ.EX length")
   call assert_true(maxval(abs(signals%blocks)) > 1.0_dp, "DJ.EX blocks signal")
   samples1 = rclaw(20, seed=8675309)
   samples2 = rclaw(20, seed=8675309)
   call assert_true(maxval(abs(samples1 - samples2)) <= tiny(1.0_dp), "seeded rclaw reproducibility")
   spectrum = cns(32)
   do level = 0, spectrum%nlevels - 1
      spectrum%detail(level)%values = 0.25_dp + 0.05_dp * real(level, dp)
   end do
   lsw_one = lsw_sim(spectrum, seed=13579)
   lsw_two = lsw_sim(spectrum, seed=13579)
   call assert_true(size(lsw_one) == 32, "LSWsim output size")
   call assert_true(maxval(abs(lsw_one - lsw_two)) <= tiny(1.0_dp), "LSWsim seeded reproducibility")
   call assert_true(maxval(abs(lsw_one)) > 0.0_dp, "LSWsim nonzero realization")

   checked_one = check_my_ews(spectrum, nsim=2, seed=24680)
   checked_two = check_my_ews(spectrum, nsim=2, seed=24680)
   call assert_true(checked_one%ok .and. checked_two%ok, "checkmyews status")
   do level = 0, checked_one%nlevels - 1
      call assert_true(maxval(abs(checked_one%detail(level)%values - checked_two%detail(level)%values)) <= &
         tiny(1.0_dp), "checkmyews seeded reproducibility")
   end do

   print *, "test_cv_signals: PASS"
contains
   subroutine assert_true(condition, message)
      logical, intent(in) :: condition !! Condition that must hold.
      character(len=*), intent(in) :: message !! Diagnostic label for a failed assertion.
      if (.not. condition) then
         print *, "FAIL: ", trim(message)
         error stop 1
      end if
   end subroutine assert_true
end program test_cv_signals
