program test_core
   use gamlss_dist
   implicit none
   real(dp) :: total, mean_value, pvalue
   integer :: k

   call assert_close(dGA(2.3_dp, 1.7_dp, 0.6_dp), &
      0.2437481154774664_dp, 2.0e-13_dp, 'dGA')
   call assert_close(pGA(2.3_dp, 1.7_dp, 0.6_dp), &
      0.7663339182045005_dp, 2.0e-13_dp, 'pGA')
   call assert_close(dWEI2(1.4_dp, 2.2_dp, 1.3_dp), &
      0.10483474557556417_dp, 2.0e-13_dp, 'dWEI2')
   call assert_close(pWEI2(1.4_dp, 2.2_dp, 1.3_dp), &
      0.9668639450531986_dp, 2.0e-13_dp, 'pWEI2')
   call assert_close(dBE(0.37_dp, 0.42_dp, 0.28_dp), &
      2.6628419108761445_dp, 2.0e-12_dp, 'dBE')
   call assert_close(pBE(0.37_dp, 0.42_dp, 0.28_dp), &
      0.3772130332749128_dp, 2.0e-13_dp, 'pBE')

   call assert_close(dNBI(4.0_dp, 2.7_dp, 0.45_dp), &
      0.10172815761673207_dp, 2.0e-13_dp, 'dNBI')
   call assert_close(pNBI(4.0_dp, 2.7_dp, 0.45_dp), &
      0.8067918215320627_dp, 2.0e-13_dp, 'pNBI')
   call assert_close(dNBII(4.0_dp, 2.7_dp, 0.45_dp), &
      0.1257594350944922_dp, 2.0e-13_dp, 'dNBII')
   call assert_close(pNBII(4.0_dp, 2.7_dp, 0.45_dp), &
      0.8314243547849951_dp, 2.0e-13_dp, 'pNBII')
   call assert_close(dPO(1.5_dp, 2.0_dp), 0.0_dp, 0.0_dp, 'dPO noninteger')
   call assert_close(dNBI(1.5_dp, 2.0_dp, 0.3_dp), 0.0_dp, 0.0_dp, 'dNBI noninteger')

   call assert_close(dBCCG(4.2_dp, 3.7_dp, 0.28_dp, 0.8_dp), &
      0.3351605223794052_dp, 3.0e-13_dp, 'dBCCG')
   call assert_close(pBCCG(4.2_dp, 3.7_dp, 0.28_dp, 0.8_dp), &
      0.6831153794299946_dp, 3.0e-13_dp, 'pBCCG')
   call assert_close(qBCCG(0.73_dp, 3.7_dp, 0.28_dp, 0.8_dp), &
      4.345418471385803_dp, 3.0e-12_dp, 'qBCCG')

   call assert_close(dBCT(4.2_dp, 3.7_dp, 0.28_dp, 0.8_dp, 5.0_dp), &
      0.31373203950875134_dp, 3.0e-13_dp, 'dBCT')
   call assert_close(pBCT(4.2_dp, 3.7_dp, 0.28_dp, 0.8_dp, 5.0_dp), &
      0.671985861893415_dp, 3.0e-13_dp, 'pBCT')
   call assert_close(qBCT(0.73_dp, 3.7_dp, 0.28_dp, 0.8_dp, 5.0_dp), &
      4.39636511467236_dp, 5.0e-12_dp, 'qBCT')

   call assert_close(dBCPE(4.2_dp, 3.7_dp, 0.28_dp, 0.8_dp, 1.5_dp), &
      0.3447276689127587_dp, 3.0e-13_dp, 'dBCPE')
   call assert_close(pBCPE(4.2_dp, 3.7_dp, 0.28_dp, 0.8_dp, 1.5_dp), &
      0.7048017244974077_dp, 3.0e-13_dp, 'pBCPE')
   call assert_close(qBCPE(0.73_dp, 3.7_dp, 0.28_dp, 0.8_dp, 1.5_dp), &
      4.275440678680146_dp, 5.0e-12_dp, 'qBCPE')

   call assert_close(pBCCG(4.2_dp, 3.7_dp, 0.28_dp, 0.0_dp), &
      0.6746120875261494_dp, 3.0e-13_dp, 'pBCCG nu=0')
   call assert_close(pBCT(4.2_dp, 3.7_dp, 0.28_dp, 0.0_dp, 5.0_dp), &
      0.6651294841662526_dp, 3.0e-13_dp, 'pBCT nu=0')
   call assert_close(pBCPE(4.2_dp, 3.7_dp, 0.28_dp, 0.0_dp, 1.5_dp), &
      0.6960479578255885_dp, 3.0e-13_dp, 'pBCPE nu=0')

   call assert_close(pBCCG(qBCCG(0.37_dp, 2.5_dp, 0.4_dp, -0.7_dp), &
      2.5_dp, 0.4_dp, -0.7_dp), 0.37_dp, 2.0e-11_dp, 'BCCG roundtrip')
   call assert_close(pBCT(qBCT(0.37_dp, 2.5_dp, 0.4_dp, -0.7_dp, 7.0_dp), &
      2.5_dp, 0.4_dp, -0.7_dp, 7.0_dp), 0.37_dp, 2.0e-10_dp, 'BCT roundtrip')
   call assert_close(pBCPE(qBCPE(0.37_dp, 2.5_dp, 0.4_dp, -0.7_dp, 1.2_dp), &
      2.5_dp, 0.4_dp, -0.7_dp, 1.2_dp), 0.37_dp, 2.0e-10_dp, 'BCPE roundtrip')

   total = 0.0_dp
   mean_value = 0.0_dp
   do k = 0, 150
      pvalue = dPIG(real(k, dp), 3.2_dp, 0.35_dp)
      total = total + pvalue
      mean_value = mean_value + real(k, dp) * pvalue
   end do
   call assert_close(total, 1.0_dp, 2.0e-11_dp, 'PIG normalization')
   call assert_close(mean_value, 3.2_dp, 3.0e-10_dp, 'PIG mean')

   total = 0.0_dp
   do k = 0, 12
      total = total + dBB(real(k, dp), 0.37_dp, 0.21_dp, 12)
   end do
   call assert_close(total, 1.0_dp, 2.0e-12_dp, 'BB normalization')

   call assert_close(pNBI(real(qNBI(0.8_dp, 2.7_dp, 0.45_dp), dp), &
      2.7_dp, 0.45_dp), 0.8_dp, 0.15_dp, 'NBI quantile coverage')

   print '(a)', 'test_core: PASS'
contains
   subroutine assert_close(actual, expected, tolerance, name)
      real(dp), intent(in) :: actual, expected, tolerance
      character(*), intent(in) :: name
      if (abs(actual - expected) > tolerance) then
         print '(a)', 'FAIL: '//trim(name)
         print '(a,es24.16)', ' actual   = ', actual
         print '(a,es24.16)', ' expected = ', expected
         error stop 1
      end if
   end subroutine assert_close
end program test_core
