! SPDX-License-Identifier: GPL-2.0-only
! Based on FinancialMath 0.1.1, Copyright (C) 2016 Kameron Penn and Jack Schmidt.
program test_annuities
   use financialmath
   use test_support
   implicit none
   type(annuity_result_t) :: a
   real(dp) :: pv, fv, direct
   integer :: j

   pv = annuity_level_pv(100.0_dp, 10.0_dp, 0.05_dp, .true.)
   call assert_close(pv, 772.1734929184818_dp)
   fv = annuity_level_fv(100.0_dp, 10.0_dp, 0.05_dp, .true.)
   call assert_close(fv, 1257.789253554884_dp)

   a = annuity_level(pv, 0.0_dp, 10.0_dp, 0.0_dp, 0.05_dp, 1.0_dp, 1.0_dp, .true., 'payment')
   call assert_true(a%status%ok, 'level annuity solve failed')
   call assert_close(a%first_payment, 100.0_dp)
   a = annuity_level(pv, 0.0_dp, 10.0_dp, 100.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, .true., 'rate')
   call assert_true(a%status%ok, 'level annuity rate solve failed')
   call assert_close(a%effective_rate, 0.05_dp, 1.0e-9_dp, 1.0e-9_dp)

   direct = 0.0_dp
   do j = 1, 8
      direct = direct+(100.0_dp+10.0_dp*real(j-1,dp))/1.04_dp**j
   end do
   call assert_close(annuity_arith_pv(100.0_dp, 10.0_dp, 8.0_dp, 0.04_dp, .true.), direct)
   a = annuity_arith(direct, 0.0_dp, 8.0_dp, 100.0_dp, 0.0_dp, 0.04_dp, 1.0_dp, 1.0_dp, .true., 'increment')
   call assert_true(a%status%ok, 'arithmetic annuity solve failed')
   call assert_close(a%increment, 10.0_dp)

   direct = 0.0_dp
   do j = 1, 12
      direct = direct+100.0_dp*1.02_dp**real(j-1,dp)/1.05_dp**j
   end do
   call assert_close(annuity_geo_pv(100.0_dp, 0.02_dp, 12.0_dp, 0.05_dp, .true.), direct)
   a = annuity_geo(direct, 0.0_dp, 12.0_dp, 100.0_dp, 0.0_dp, 0.05_dp, 1.0_dp, 1.0_dp, .true., 'growth_rate')
   call assert_true(a%status%ok, 'geometric annuity growth solve failed')
   call assert_close(a%growth_rate, 0.02_dp, 1.0e-9_dp, 1.0e-9_dp)

   a = perpetuity_level(0.0_dp, 100.0_dp, 0.05_dp, 1.0_dp, 1.0_dp, .true., 'pv')
   call assert_true(a%status%ok, 'level perpetuity failed')
   call assert_close(a%present_value, 2000.0_dp)
   a = perpetuity_arith(0.0_dp, 100.0_dp, 5.0_dp, 0.05_dp, 1.0_dp, 1.0_dp, .true., 'pv')
   call assert_close(a%present_value, 4000.0_dp)
   a = perpetuity_geo(0.0_dp, 100.0_dp, 0.02_dp, 0.05_dp, 1.0_dp, 1.0_dp, .true., 'pv')
   call assert_close(a%present_value, 100.0_dp/0.03_dp)

   print '(a)', 'test_annuities: PASS'
end program test_annuities
