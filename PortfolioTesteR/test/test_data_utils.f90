! SPDX-License-Identifier: MIT
! PortfolioTesteR modern Fortran translation
program test_data_utils
  use portfolio_tester
  implicit none
  real(dp)::a(4,2)
  real(dp),allocatable::ff(:,:),aligned(:,:),z(:,:)
  integer::src(4),dst(3)
  a=reshape([1.0_dp,nan_dp(),3.0_dp,4.0_dp,10.0_dp,11.0_dp,nan_dp(),13.0_dp],[4,2])
  call forward_fill(a,ff)
  call assert_close(ff(2,1),1.0_dp,1.0e-12_dp,'forward fill')
  src=[1,3,5,7];dst=[2,5,8]
  call align_to_indices(src,a,dst,aligned,.true.)
  call assert_close(aligned(2,1),a(3,1),1.0e-12_dp,'alignment exact')
  call standardize_panel(ff,z,.false.)
  call assert_true(count(is_finite(z))>0,'standardization')
  print '(a)','test_data_utils: PASS'
contains
  subroutine assert_true(ok,msg)
    logical,intent(in)::ok;character(len=*),intent(in)::msg
    if(.not.ok)then;write(*,'(a)')'FAIL: '//msg;error stop 1;end if
  end subroutine
  subroutine assert_close(x,y,tol,msg)
    real(dp),intent(in)::x,y,tol;character(len=*),intent(in)::msg
    call assert_true(abs(x-y)<=tol,msg)
  end subroutine
end program test_data_utils
