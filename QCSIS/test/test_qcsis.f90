program test_qcsis
   use qcsis_mod, only : cqc, cqcsis, dp, qc, qc_result, qcsis, screening_result
   implicit none

   real(dp), parameter :: tol = 5.0e-13_dp
   real(dp) :: x1(6), y(6), xmat(6, 3), xtie(6, 2), cqc_value
   real(dp) :: tau(3)
   type(qc_result) :: qfit
   type(screening_result) :: sfit
   integer :: stat
   character(len=:), allocatable :: errmsg

   x1 = [1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, 16.0_dp, 32.0_dp]
   y = [3.0_dp, 1.0_dp, 4.0_dp, 1.5_dp, 5.0_dp, 9.0_dp]
   tau = [0.25_dp, 0.50_dp, 0.75_dp]

   xmat(:, 1) = x1
   xmat(:, 2) = [2.0_dp, 1.0_dp, 3.0_dp, 6.0_dp, 5.0_dp, 4.0_dp]
   xmat(:, 3) = [6.0_dp, 5.0_dp, 4.0_dp, 3.0_dp, 2.0_dp, 1.0_dp]

   qfit = qc(x1, y, tau, stat, errmsg)
   call assert_true(stat == 0, "qc status")
   call assert_close(qfit%rho(1), 0.356938792420570_dp, tol, "qc rho(1)")
   call assert_close(qfit%rho(2), 0.576083660687569_dp, tol, "qc rho(2)")
   call assert_close(qfit%rho(3), 0.876122490486854_dp, tol, "qc rho(3)")

   cqc_value = cqc(x1, y, stat, errmsg)
   call assert_true(stat == 0, "cqc status")
   call assert_close(cqc_value, 0.567958408816597_dp, tol, "cqc value")

   sfit = qcsis(xmat, y, tau, 2, stat, errmsg)
   call assert_true(stat == 0, "qcsis explicit tau status")
   call assert_close(sfit%w(1), 0.408956101327577_dp, tol, "qcsis w(1)")
   call assert_close(sfit%w(2), 0.080246913580247_dp, tol, "qcsis w(2)")
   call assert_close(sfit%w(3), 0.369488536155203_dp, tol, "qcsis w(3)")
   call assert_true(all(sfit%selected == [1, 3]), "qcsis selected")

   sfit = cqcsis(xmat, y, 2, stat, errmsg)
   call assert_true(stat == 0, "cqcsis status")
   call assert_close(sfit%w(1), 0.567958408816597_dp, tol, "cqcsis w(1)")
   call assert_close(sfit%w(2), 0.272472576104458_dp, tol, "cqcsis w(2)")
   call assert_close(sfit%w(3), 0.504940726890581_dp, tol, "cqcsis w(3)")
   call assert_true(all(sfit%selected == [1, 3]), "cqcsis selected")

   sfit = qcsis(xmat, y, 2, stat, errmsg)
   call assert_true(stat == 0, "qcsis default tau status")

   xtie(:, 1) = x1
   xtie(:, 2) = x1
   sfit = qcsis(xtie, y, tau, 2, stat, errmsg)
   call assert_true(stat == 0, "qcsis tied predictors status")
   call assert_true(all(sfit%selected == [1, 2]), "ties preserve predictor order")

   qfit = qc(x1, y, [0.0_dp], stat, errmsg)
   call assert_true(stat /= 0, "invalid tau rejected")
   call assert_true(size(qfit%rho) == 0, "failed qc has empty result")

   print '(a)', "All QCSIS tests passed."

contains

   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label

      if (abs(actual - expected) > tolerance) then
         print '(a)', "FAIL: " // label
         print '(a,es24.16)', "  actual:   ", actual
         print '(a,es24.16)', "  expected: ", expected
         error stop 1
      end if
   end subroutine assert_close


   subroutine assert_true(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label

      if (.not. condition) then
         print '(a)', "FAIL: " // label
         error stop 1
      end if
   end subroutine assert_true

end program test_qcsis
