! SPDX-License-Identifier: GPL-3.0-or-later
program test_distributions
  use degreenet_kinds, only : dp
  use degreenet_math, only : zeta_r
  use degreenet_distributions, only : dyule, ddp, dwar, ddqe, dnb_degree
  implicit none
  real(dp) :: v2(2)
  call check_close(zeta_r(2.0_dp),1.6449340668482264_dp,2e-13_dp,'zeta(2)')
  call check_close(dyule(3.5_dp,2),0.15873015873015872_dp,2e-14_dp,'Yule')
  call check_close(ddp(3.5_dp,2),0.07844651715208137_dp,2e-12_dp,'Zipf')
  v2=[3.5_dp,0.4_dp]
  call check_close(dwar(v2,2),0.18315018315018305_dp,2e-14_dp,'Waring')
  v2=[3.5_dp,1.2_dp]
  call check_close(ddqe(v2,2),0.08756159724936374_dp,2e-14_dp,'DQE')
  v2=[5.0_dp,0.2_dp]
  call check_close(dnb_degree(v2,0,0,.false.),0.2_dp,2e-14_dp,'NB p0')
  print *, 'test_distributions: PASS'
contains
  subroutine check_close(a,b,tol,name)
    real(dp),intent(in)::a,b,tol;character(*),intent(in)::name
    if(abs(a-b)>tol*max(1.0_dp,abs(b)))then
      print *, 'FAIL ',name,a,b;error stop 1
    end if
  end subroutine
end program
