! SPDX-License-Identifier: GPL-2.0-or-later
program test_continuous
  use actuar, only : dp, dpareto, ppareto, qpareto, dburr, pburr, &
    pgenpareto, pllogis, pinvgamma, pinvweibull, ptrgamma, &
    pgenbeta, pgumbel, pinvgauss
  implicit none

  call assert_close(dpareto(4.0_dp,2.0_dp,3.0_dp), &
    0.0524781341107872_dp,1.0e-13_dp)
  call assert_close(ppareto(4.0_dp,2.0_dp,3.0_dp), &
    0.8163265306122449_dp,1.0e-13_dp)
  call assert_close(qpareto(0.8_dp,2.0_dp,3.0_dp), &
    3.708203932499369_dp,1.0e-12_dp)
  call assert_close(dburr(2.0_dp,1.7_dp,2.3_dp,4.0_dp), &
    0.2409900296920668_dp,1.0e-12_dp)
  call assert_close(pburr(2.0_dp,1.7_dp,2.3_dp,4.0_dp), &
    0.2696858595391708_dp,1.0e-12_dp)
  call assert_close(pgenpareto(5.0_dp,1.4_dp,2.2_dp,3.0_dp), &
    0.4874277830716066_dp,2.0e-12_dp)
  call assert_close(pllogis(2.0_dp,1.5_dp,3.0_dp), &
    0.3524704450894247_dp,1.0e-13_dp)
  call assert_close(pinvgamma(2.0_dp,3.0_dp,4.0_dp), &
    0.6766764161830634_dp,2.0e-12_dp)
  call assert_close(pinvweibull(2.0_dp,2.0_dp,3.0_dp), &
    0.1053992245618643_dp,1.0e-13_dp)
  call assert_close(ptrgamma(2.0_dp,2.0_dp,1.5_dp,3.0_dp), &
    0.1039331076481856_dp,2.0e-12_dp)
  call assert_close(pgenbeta(2.0_dp,2.0_dp,3.0_dp,1.4_dp,5.0_dp), &
    0.3084527822167539_dp,2.0e-12_dp)
  call assert_close(pgumbel(1.0_dp,0.5_dp,2.0_dp), &
    0.4589560693076638_dp,1.0e-13_dp)
  call assert_close(pinvgauss(1.5_dp,2.0_dp,0.3_dp), &
    0.4820622969411085_dp,3.0e-12_dp)

  print '(a)', 'test_continuous: PASS'
contains
  subroutine assert_close(actual,expected,tol)
    real(dp), intent(in) :: actual,expected,tol
    if(abs(actual-expected)>tol*max(1.0_dp,abs(expected))) then
      print '(a,3es24.15)', 'mismatch: ',actual,expected,abs(actual-expected)
      error stop 1
    end if
  end subroutine assert_close
end program test_continuous
