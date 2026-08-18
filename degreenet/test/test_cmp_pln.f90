! SPDX-License-Identifier: GPL-3.0-or-later
program test_cmp_pln
  use degreenet_kinds, only : dp
  use degreenet_distributions, only : dcmp_natural, cmp_moments, cmp_mutonatural, dpln
  implicit none
  real(dp)::mu,sd,lambda,nu
  logical::ok
  call check_close(dcmp_natural(2.0_dp,1.0_dp,3),0.1804470443154836_dp,2e-12_dp,'CMP=Poisson')
  call cmp_moments(2.0_dp,1.0_dp,mu,sd)
  call check_close(mu,2.0_dp,2e-9_dp,'CMP mean')
  call check_close(sd,sqrt(2.0_dp),2e-9_dp,'CMP sd')
  call cmp_mutonatural(2.0_dp,sqrt(2.0_dp),lambda,nu,ok)
  if(.not.ok)then;print *,'FAIL CMP conversion';error stop 1;end if
  call check_close(lambda,2.0_dp,2e-5_dp,'CMP lambda inversion')
  call check_close(nu,1.0_dp,2e-5_dp,'CMP nu inversion')
  call check_close(dpln([0.0_dp,0.5_dp],2),0.17717630953297012_dp,2e-8_dp,'PLN')
  print *, 'test_cmp_pln: PASS'
contains
  subroutine check_close(a,b,tol,name)
    real(dp),intent(in)::a,b,tol;character(*),intent(in)::name
    if(abs(a-b)>tol*max(1.0_dp,abs(b)))then
      print *, 'FAIL ',name,a,b;error stop 1
    end if
  end subroutine
end program
