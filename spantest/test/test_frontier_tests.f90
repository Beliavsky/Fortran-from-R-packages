! SPDX-License-Identifier: GPL-3.0-only
program test_frontier_tests
  use spantest
  implicit none
  type(simulation_result) :: sim
  type(span_result) :: f2,hk,km,bad
  real(dp) :: tiny_r1(4,3),tiny_r2(4,2)
  integer :: i

  sim=span_simulate(220,4,3,ncp=0.0_dp,dgp=4,seed=19)
  call check(sim%status==span_ok,'simulation failed')
  f2=span_f2(sim%r1,sim%r2)
  hk=span_hk(sim%r1,sim%r2)
  km=span_km(sim%r1,sim%r2)
  call check(f2%status==span_ok .and. hk%status==span_ok .and. km%status==span_ok,'frontier test failed')
  call check(in_unit(f2%pval) .and. in_unit(hk%pval) .and. in_unit(km%pval),'invalid frontier p-value')
  call check(f2%stat>-1.0e-10_dp .and. hk%stat>-1.0e-10_dp .and. km%stat>-1.0e-10_dp,'invalid statistic')

  tiny_r1=reshape([(real(i,dp),i=1,12)],[4,3])
  tiny_r2=reshape([(real(i,dp),i=1,8)],[4,2])
  bad=span_f2(tiny_r1,tiny_r2)
  call check(bad%status==span_insufficient_df,'degrees-of-freedom guard failed')
  print '(a)','test_frontier_tests: PASS'
contains
  pure logical function in_unit(x)
    real(dp),intent(in)::x
    in_unit=x>=0.0_dp .and. x<=1.0_dp
  end function
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if (.not.ok) error stop msg
  end subroutine
end program test_frontier_tests
