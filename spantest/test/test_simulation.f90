! SPDX-License-Identifier: GPL-3.0-only
program test_simulation
  use spantest
  implicit none
  type(simulation_result) :: a,b,s
  real(dp) :: z(5),e(5),expected(5),s2
  integer :: d,i

  do d=1,12
    s=span_simulate(80,3,4,ncp=0.1_dp,dgp=d,seed=100+d,burnin=50)
    call check(s%status==span_ok,'DGP preset failed')
    call check(all(shape(s%r1)==[80,3]) .and. all(shape(s%r2)==[80,4]),'wrong simulation dimensions')
  end do
  a=span_simulate(100,2,3,ncp=0.2_dp,dgp=9,seed=55,burnin=30)
  b=span_simulate(100,2,3,ncp=0.2_dp,dgp=9,seed=55,burnin=30)
  call check(maxval(abs(a%r1-b%r1))<tiny(1.0_dp) .and. &
             maxval(abs(a%r2-b%r2))<tiny(1.0_dp),'simulation is not reproducible')

  z=[1.0_dp,-0.5_dp,0.25_dp,2.0_dp,-1.0_dp]
  call garch_filter(z,0.1_dp,0.1_dp,0.8_dp,e)
  s2=1.0_dp
  do i=1,5
    expected(i)=sqrt(s2)*z(i)
    s2=0.1_dp+0.1_dp*expected(i)**2+0.8_dp*s2
  end do
  call check(maxval(abs(e-expected))<1.0e-14_dp,'GARCH recursion mismatch')
  s=span_simulate(10,2,2,dgp=13)
  call check(s%status==span_invalid_input,'invalid DGP not rejected')
  print '(a)','test_simulation: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if (.not.ok) error stop msg
  end subroutine
end program test_simulation
