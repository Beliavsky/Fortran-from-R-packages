! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Modern Fortran translation of computational code from the R package sn.
program test_rng_moments
  use iso_fortran_env, only : int64
  use sn, only : dp, sn_rng_state, rsn, rst, sn_uv_params, sn_cumulants, sn_ok
  implicit none
  type(sn_rng_state) :: rng
  type(sn_uv_params) :: par
  real(dp),allocatable :: x(:),cum(:)
  real(dp) :: mean,var
  integer :: info
  allocate(x(40000))
  call rng%seed(123456789_int64)
  call rsn(rng,x,xi=0.5_dp,omega=1.7_dp,alpha=2.0_dp,info=info)
  call assert_true(info==sn_ok,'rsn status')
  par=sn_uv_params(0.5_dp,1.7_dp,2.0_dp,0.0_dp)
  call sn_cumulants(par,2,cum,info)
  mean=sum(x)/real(size(x),dp)
  var=sum((x-mean)**2)/real(size(x)-1,dp)
  call assert_close(mean,cum(1),0.025_dp,'sample mean')
  call assert_close(var,cum(2),0.04_dp,'sample variance')
  call rng%seed(998877_int64)
  call rst(rng,x,alpha=-1.0_dp,nu=6.0_dp,info=info)
  call assert_true(info==sn_ok .and. all(abs(x)<1.0e6_dp),'rst finite')
  print '(a)','test_rng_moments: PASS'
contains
  subroutine assert_close(a,b,tol,msg)
    real(dp),intent(in)::a,b,tol
    character(len=*),intent(in)::msg
    if(abs(a-b)>tol*max(1.0_dp,abs(a),abs(b))) then
      write(*,*) trim(msg),a,b
      error stop 1
    end if
  end subroutine
  subroutine assert_true(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) then
      write(*,*) trim(msg)
      error stop 1
    end if
  end subroutine
end program test_rng_moments
