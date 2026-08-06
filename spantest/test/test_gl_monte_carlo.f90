! SPDX-License-Identifier: GPL-3.0-only
program test_gl_monte_carlo
  use spantest
  implicit none
  type(simulation_result) :: sim0,sim1
  type(gl_result) :: a1,a2,j1

  sim0=span_simulate(220,3,5,ncp=0.0_dp,dgp=7,seed=7)
  sim1=span_simulate(220,3,5,ncp=0.5_dp,dgp=7,seed=7)
  call check(sim0%status==span_ok .and. sim1%status==span_ok,'simulation failed')
  a1=span_gl_a(sim0%r1,sim0%r2,totsim=79,seed=91)
  a2=span_gl_a(sim0%r1,sim0%r2,totsim=79,seed=91)
  call check(a1%status==span_ok .and. a2%status==span_ok,'GL alpha test failed')
  call check(abs(a1%pval_lmc-a2%pval_lmc)<tiny(1.0_dp) .and. &
             abs(a1%pval_bmc-a2%pval_bmc)<tiny(1.0_dp),'GL seed is not reproducible')
  call check(in_unit(a1%pval_lmc) .and. in_unit(a1%pval_bmc),'invalid GL p-value')
  call check(valid_decision(a1%decision),'invalid GL decision')

  j1=span_gl_ad(sim1%r1,sim1%r2,totsim=99,seed=3)
  call check(j1%status==span_ok,'GL joint test failed')
  call check(j1%pval_bmc<=0.05_dp,'strong joint alternative not rejected')
  print '(a)','test_gl_monte_carlo: PASS'
contains
  pure logical function in_unit(x)
    real(dp),intent(in)::x
    in_unit=x>=0.0_dp .and. x<=1.0_dp
  end function
  pure logical function valid_decision(i)
    integer,intent(in)::i
    valid_decision=i==-1 .or. i==0 .or. i==1
  end function
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if (.not.ok) error stop msg
  end subroutine
end program test_gl_monte_carlo
