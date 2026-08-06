! SPDX-License-Identifier: GPL-3.0-only
program test_subseries_cauchy
  use spantest
  implicit none
  type(simulation_result) :: sim0,sim1
  type(as_result) :: a0,a0b,a1
  integer :: i

  sim0=span_simulate(250,3,6,ncp=0.0_dp,dgp=1,seed=3)
  sim1=span_simulate(250,3,6,ncp=0.5_dp,dgp=1,seed=3)
  a0=span_as(sim0%r1,sim0%r2,b_draws=5,seed=22)
  a0b=span_as(sim0%r1,sim0%r2,b_draws=5,seed=22)
  call check(a0%status==span_ok .and. a0b%status==span_ok,'subseries test failed')
  call check(size(a0%pvalues)==6,'unexpected default result count')
  call check(maxval(abs(a0%pvalues-a0b%pvalues))<tiny(1.0_dp),'subseries seed is not reproducible')
  do i=1,size(a0%pvalues)
    call check(a0%pvalues(i)>=0.0_dp .and. a0%pvalues(i)<=1.0_dp,'invalid Cauchy p-value')
  end do
  call check(trim(a0%names(1))=='CCTd_L0_k1','unexpected first result name')
  call check(trim(a0%names(6))=='CCTa_L2_k1','unexpected last result name')

  a1=span_as(sim1%r1,sim1%r2,ks=[1.0_dp/3.0_dp],l_values=[0],seed=3)
  call check(a1%status==span_ok,'alternative subseries test failed')
  call check(trim(a1%names(3))=='CCTa_L0_k1','unexpected alpha result position')
  call check(a1%pvalues(3)<0.01_dp,'alpha alternative not detected')
  print '(a)','test_subseries_cauchy: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if (.not.ok) error stop msg
  end subroutine
end program test_subseries_cauchy
