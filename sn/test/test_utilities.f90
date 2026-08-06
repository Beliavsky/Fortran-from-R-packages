! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Modern Fortran translation of computational code from the R package sn.
program test_utilities
  use sn, only : dp, pi, vech, vech_to_matrix, duplication_matrix, block_diag, trace, &
                 pprodn2, pprodt2, qprodt2, q_penalty, sn_ok, fit_grouped, grouped_fit_result
  implicit none
  real(dp) :: a(2,2),b(1,1),expected,p,q
  real(dp),allocatable :: v(:),aa(:,:),d(:,:),c(:,:)
  real(dp) :: lo(7),hi(7),cnt(7)
  type(grouped_fit_result) :: fit
  integer :: info
  a=reshape([2.0_dp,0.5_dp,0.5_dp,3.0_dp],[2,2])
  call vech(a,v,info)
  call vech_to_matrix(v,aa,info)
  call assert_true(info==sn_ok .and. maxval(abs(aa-a))<1.0e-14_dp,'vech')
  call duplication_matrix(2,d,info)
  call assert_true(size(d,1)==4 .and. size(d,2)==3,'duplication size')
  b(1,1)=4.0_dp
  call block_diag(a,b,c)
  call assert_close(trace(c),9.0_dp,1.0e-14_dp,'block trace')
  expected=0.5_dp-asin(0.35_dp)/pi
  call assert_close(pprodn2(0.0_dp,0.35_dp),expected,2.0e-5_dp,'normal product sign')
  call assert_close(pprodt2(0.0_dp,0.35_dp,7.0_dp),expected,3.0e-3_dp,'t product sign')
  p=0.62_dp
  q=qprodt2(p,-0.2_dp,8.0_dp,info=info)
  call assert_true(info==sn_ok,'qprodt status')
  call assert_close(pprodt2(q,-0.2_dp,8.0_dp),p,2.0e-4_dp,'qprodt inversion')
  call assert_true(q_penalty(2.0_dp)>0.0_dp,'Q penalty')

  lo=[-3.5_dp,-2.5_dp,-1.5_dp,-0.5_dp,0.5_dp,1.5_dp,2.5_dp]
  hi=lo+1.0_dp
  cnt=[1.0_dp,6.0_dp,24.0_dp,38.0_dp,24.0_dp,6.0_dp,1.0_dp]
  call fit_grouped(lo,hi,cnt,'NORMAL',fit,max_iter=1800,tol=1.0e-8_dp)
  call assert_true(fit%status==sn_ok .and. fit%omega>0.7_dp .and. fit%omega<1.5_dp,'grouped fit')
  print '(a)','test_utilities: PASS'
contains
  subroutine assert_close(x,y,tol,msg)
    real(dp),intent(in)::x,y,tol
    character(len=*),intent(in)::msg
    if(abs(x-y)>tol*max(1.0_dp,abs(x),abs(y))) then
      write(*,*) trim(msg),x,y
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
end program test_utilities
