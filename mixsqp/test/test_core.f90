program test_core
  use mixsqp
  implicit none
  integer, parameter :: n=6,m=4
  real(dp) :: L(n,m),w(n),xref(m)
  type(mixsqp_result) :: r
  type(mixsqp_control) :: c
  L = reshape([ &
    1.0_dp,0.2_dp,0.8_dp,0.3_dp,1.1_dp,0.5_dp, &
    0.3_dp,1.2_dp,0.5_dp,0.7_dp,0.2_dp,0.9_dp, &
    0.6_dp,0.4_dp,1.3_dp,0.2_dp,0.8_dp,0.3_dp, &
    0.2_dp,0.7_dp,0.4_dp,1.1_dp,0.5_dp,1.0_dp], [n,m])
  w=[1._dp,2._dp,1._dp,1._dp,3._dp,2._dp]
  xref=[0.485912573_dp,0.195013132_dp,0._dp,0.319074295_dp]
  c=mixsqp_default_control(); c%tol_svd=0._dp; c%verbose=.false.
  call fit_mixsqp(L,r,w_in=w,control=c)
  call check(r%status==0,'solver did not converge')
  call check(maxval(abs(r%x-xref))<5e-7_dp,'solution mismatch')
  call check(abs(r%value-0.4161754053877256_dp)<2e-10_dp,'objective mismatch')
  call check(abs(sum(r%x)-1._dp)<1e-12_dp,'solution not normalized')
  call check(minval(r%x)>=0._dp,'solution infeasible')
  print *, 'test_core: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) error stop msg
  end subroutine
end program
