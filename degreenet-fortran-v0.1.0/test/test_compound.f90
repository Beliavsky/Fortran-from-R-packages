! SPDX-License-Identifier: GPL-3.0-or-later
program test_compound
  use degreenet_kinds, only : dp
  use degreenet_compound, only : dgyule, dnbyule, dgwar, waring_prob_to_natural, nbmean
  implicit none
  real(dp)::a,gm,gs,s
  integer::k
  call check_close(dgyule([3.5_dp,5.0_dp],2),0.11784126984126986_dp,2e-13_dp,'geometric-Yule')
  call waring_prob_to_natural(3.5_dp,0.2_dp,a)
  call check_close(a,5.0_dp,1e-14_dp,'Waring conversion')
  call nbmean([5.0_dp,0.2_dp],gm,gs)
  call check_close(gm,4.0_dp,1e-14_dp,'NB gamma mean')
  call check_close(gs,4.0_dp,1e-14_dp,'NB gamma sd')
  s=0.0_dp;do k=1,2000;s=s+dgyule([3.5_dp,5.0_dp],k);end do
  call check_close(s,0.8_dp,5e-7_dp,'geometric-Yule positive mass')
  if(dnbyule([3.5_dp,5.0_dp,0.2_dp],2)<=0.0_dp)error stop 1
  if(dgwar([3.5_dp,0.4_dp,5.0_dp],2)<=0.0_dp)error stop 1
  print *, 'test_compound: PASS'
contains
  subroutine check_close(x,y,tol,name)
    real(dp),intent(in)::x,y,tol;character(*),intent(in)::name
    if(abs(x-y)>tol*max(1.0_dp,abs(y)))then;print *,'FAIL ',name,x,y;error stop 1;end if
  end subroutine
end program
