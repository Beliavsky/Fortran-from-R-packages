program test_core
   use risk_parity_portfolio_mod
   implicit none
   integer, parameter :: n = 8
   real(dp) :: vol(n), corr(n, n), sigma(n, n), b(n)
   real(dp) :: w_spinu(n), w_roncalli(n), w_choi(n), w_newton(n), w_active(n)
   real(dp) :: expected_w(n), mu(n), diag_sigma(3, 3), diag_b(3), diag_w(3)
   integer :: i, j, info, iterations, failures
   type(risk_parity_result) :: result

   failures = 0
   vol = [0.05_dp, 0.05_dp, 0.07_dp, 0.10_dp, 0.15_dp, 0.15_dp, 0.15_dp, 0.18_dp]
   corr = reshape([ &
      100, 80, 60,-20,-10,-20,-20,-20, &
       80,100, 40,-20,-20,-10,-20,-20, &
       60, 40,100, 50, 30, 20, 20, 30, &
      -20,-20, 50,100, 60, 60, 50, 60, &
      -10,-20, 30, 60,100, 90, 70, 70, &
      -20,-10, 20, 60, 90,100, 60, 70, &
      -20,-20, 20, 50, 70, 60,100, 70, &
      -20,-20, 30, 60, 70, 70, 70,100], [n, n]) / 100.0_dp
   do j = 1, n
      do i = 1, n
         sigma(i, j) = corr(i, j) * vol(i) * vol(j)
      end do
   end do
   b = 1.0_dp / real(n, dp)
   expected_w = [0.268306_dp, 0.286769_dp, 0.114095_dp, 0.097985_dp, &
                 0.056135_dp, 0.059029_dp, 0.066560_dp, 0.051121_dp]

   call risk_parity_ccd_spinu(sigma, b, w_spinu, info, iterations, 1.0e-11_dp, 10000)
   call assert_true(info == RPP_OK, 'Spinu solver status', failures)
   call assert_close(maxval(abs(w_spinu - expected_w)), 0.0_dp, 2.0e-6_dp, &
                     'published pyrb weights', failures)
   call assert_close(maxval(abs(relative_risk_contributions(sigma, w_spinu) - b)), &
                     0.0_dp, 1.0e-9_dp, 'Spinu risk budgets', failures)

   call risk_parity_ccd_roncalli(sigma, b, w_roncalli, info, iterations, 1.0e-11_dp, 10000)
   call assert_true(info == RPP_OK, 'Roncalli solver status', failures)
   call risk_parity_ccd_choi(sigma, b, w_choi, info, iterations, 1.0e-11_dp, 10000)
   call assert_true(info == RPP_OK, 'Choi solver status', failures)
   call risk_parity_newton(sigma, b, w_newton, info, iterations, 1.0e-11_dp, 1000)
   call assert_true(info == RPP_OK, 'Newton solver status', failures)
   call assert_close(maxval(abs(w_roncalli - w_spinu)), 0.0_dp, 2.0e-8_dp, &
                     'Roncalli versus Spinu', failures)
   call assert_close(maxval(abs(w_choi - w_spinu)), 0.0_dp, 2.0e-8_dp, &
                     'Choi versus Spinu', failures)
   call assert_close(maxval(abs(w_newton - w_spinu)), 0.0_dp, 2.0e-8_dp, &
                     'Newton versus Spinu', failures)

   mu = [0.04_dp, 0.03_dp, 0.05_dp, 0.06_dp, 0.07_dp, 0.065_dp, 0.055_dp, 0.045_dp]
   call active_risk_parity_ccd(sigma, b, mu, 1.0e6_dp, 0.0_dp, w_active, info, &
                               iterations, 1.0e-6_dp, 10000)
   call assert_true(info == RPP_OK, 'active solver status', failures)
   call assert_close(maxval(abs(w_active - w_roncalli)), 0.0_dp, 2.0e-4_dp, &
                     'active solver large-tradeoff limit', failures)

   diag_sigma = 0.0_dp
   diag_sigma(1, 1) = 0.04_dp
   diag_sigma(2, 2) = 0.09_dp
   diag_sigma(3, 3) = 0.16_dp
   diag_b = [0.2_dp, 0.3_dp, 0.5_dp]
   diag_w = diagonal_risk_parity(diag_sigma, diag_b)
   call assert_close(sum(diag_w), 1.0_dp, 1.0e-14_dp, 'diagonal weights budget', failures)
   call assert_close(maxval(abs(relative_risk_contributions(diag_sigma, diag_w) - diag_b)), &
                     0.0_dp, 1.0e-13_dp, 'diagonal risk budgets', failures)

   call risk_parity_portfolio(sigma, result, b=b)
   call assert_true(result%status == RPP_OK, 'high-level status', failures)
   call assert_true(result%feasible, 'high-level feasibility', failures)
   call assert_close(maxval(abs(result%weights - expected_w)), 0.0_dp, 2.0e-6_dp, &
                     'high-level weights', failures)
   call assert_close(sum(result%weights), 1.0_dp, 1.0e-13_dp, 'high-level budget', failures)

   if (failures > 0) then
      write(*, '(a,i0)') 'test_core failures: ', failures
      error stop 1
   end if
   write(*, '(a)') 'test_core: all tests passed'
contains
   subroutine assert_true(condition, name, failures)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: name
      integer, intent(inout) :: failures
      if (.not. condition) then
         failures = failures + 1
         write(*, '(a)') 'FAIL: ' // trim(name)
      end if
   end subroutine assert_true

   subroutine assert_close(value, target, tolerance, name, failures)
      real(dp), intent(in) :: value, target, tolerance
      character(len=*), intent(in) :: name
      integer, intent(inout) :: failures
      call assert_true(abs(value - target) <= tolerance, name, failures)
   end subroutine assert_close
end program test_core
