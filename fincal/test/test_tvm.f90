! SPDX-License-Identifier: GPL-2.0-or-later
program test_tvm
   use fincal
   implicit none
   real(dp), parameter :: tol = 2.0e-10_dp
   real(dp), parameter :: cf(4) = [-5.0_dp, 1.6_dp, 2.4_dp, 2.8_dp]
   real(dp) :: rate
   integer :: status

   call assert_close(fv_simple(0.08_dp, 10.0_dp, -300.0_dp), 647.677499181836_dp, 1.0e-9_dp, 'fv_simple')
   call assert_close(pv_simple(0.07_dp, 10.0_dp, 100.0_dp), -50.834929213472_dp, 1.0e-6_dp, 'pv_simple')
   call assert_close(fv_annuity(0.03_dp, 12.0_dp, -1000.0_dp), 14192.02956_dp, 1.0e-5_dp, 'fv_annuity')
   call assert_close(pv_annuity(0.03_dp, 12.0_dp, 1000.0_dp), -9954.00399356757_dp, 1.0e-5_dp, 'pv_annuity')
   call assert_close(fv(0.07_dp, 10.0_dp, present_value = 1000.0_dp, payment = 10.0_dp), &
      -2105.31583690236_dp, 1.0e-9_dp, 'fv')
   call assert_close(pv(0.07_dp, 10.0_dp, future_value = 1000.0_dp, payment = 10.0_dp), &
      -578.585107544044_dp, 1.0e-9_dp, 'pv')
   call assert_close(fv_uneven(0.1_dp, [-1000.0_dp, -500.0_dp, 0.0_dp, 4000.0_dp, 3500.0_dp, 2000.0_dp]), &
      -8347.44_dp, 1.0e-9_dp, 'fv_uneven')
   call assert_close(pv_uneven(0.1_dp, [-1000.0_dp, -500.0_dp, 0.0_dp, 4000.0_dp, 3500.0_dp, 2000.0_dp]), &
      -4711.91226268810_dp, 1.0e-9_dp, 'pv_uneven')
   call assert_close(npv(0.12_dp, cf), 0.334821428571428_dp, tol, 'npv')

   rate = irr(cf, status)
   call assert_equal_int(status, fincal_ok, 'irr status')
   call assert_close(npv(rate, cf), 0.0_dp, 1.0e-10_dp, 'irr root')
   call assert_close(rate, 0.155175727575292_dp, 1.0e-9_dp, 'irr value')

   rate = irr2(cf, cutoff = 1.0e-3_dp, from_rate = 0.15_dp, to_rate = 0.16_dp, &
      step = 1.0e-4_dp, status = status)
   call assert_equal_int(status, fincal_ok, 'irr2 status')
   call assert_close(rate, 0.1551_dp, 1.0e-12_dp, 'irr2 value')

   rate = discount_rate(5.0_dp, 0.0_dp, 600.0_dp, -100.0_dp, status = status)
   call assert_equal_int(status, fincal_ok, 'discount_rate status')
   call assert_close(rate, 0.0912806233_dp, 1.0e-9_dp, 'discount_rate')
   call assert_close(pmt(0.08_dp, 10.0_dp, -1000.0_dp, 0.0_dp), 149.0294887_dp, 1.0e-7_dp, 'pmt')
   call assert_close(n_period(0.1_dp, -10000.0_dp, 60000000.0_dp, -50000.0_dp), 50.1099455284402_dp, 1.0e-7_dp, 'n_period')
   call assert_close(pv_perpetuity(0.1_dp, 1000.0_dp, 0.02_dp), -12500.0_dp, tol, 'pv_perpetuity')
   call assert_close(perpetuity_rate(4.5_dp, -75.0_dp), 0.06_dp, tol, 'perpetuity_rate')
   call assert_close(r_perpetuity(4.5_dp, -75.0_dp), 0.06_dp, tol, 'r_perpetuity alias')

   rate = irr([1.0_dp, 2.0_dp, 3.0_dp], status)
   call assert_equal_int(status, fincal_no_root, 'irr no-root status')
   print '(a)', 'test_tvm: PASS'
contains
   subroutine assert_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: message
      if (abs(actual - expected) > tolerance) then
         print *, trim(message), actual, expected, abs(actual - expected)
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_equal_int(actual, expected, message)
      integer, intent(in) :: actual, expected
      character(len=*), intent(in) :: message
      if (actual /= expected) then
         print *, trim(message), actual, expected
         error stop 1
      end if
   end subroutine assert_equal_int
end program test_tvm
