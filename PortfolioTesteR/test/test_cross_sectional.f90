! SPDX-License-Identifier: MIT
! PortfolioTesteR modern Fortran translation
program test_cross_sectional
  use portfolio_tester
  implicit none
  real(dp)::x(2,4),cond(2,4)
  integer::groups(4)
  real(dp),allocatable::r(:,:),breadth(:),gb(:,:),rel(:,:)
  x(1,:)=[1.0_dp,4.0_dp,2.0_dp,3.0_dp]
  x(2,:)=[2.0_dp,1.0_dp,4.0_dp,3.0_dp]
  cond=merge(1.0_dp,0.0_dp,x>2.0_dp);groups=[1,1,2,2]
  call calc_relative_strength_rank(x,r,'percentile')
  call assert_true(all(r>=0.0_dp.and.r<=1.0_dp),'rank range')
  call calc_market_breadth(cond,breadth)
  call assert_close(breadth(1),0.5_dp,1.0e-12_dp,'breadth row1')
  call group_breadth(cond,groups,gb)
  call assert_true(all(gb>=0.0_dp.and.gb<=1.0_dp),'group breadth')
  call group_relative_indicators(x,groups,rel)
  call assert_close(rel(1,1)+rel(1,2),0.0_dp,1.0e-12_dp,'group centered 1')
  call assert_close(rel(1,3)+rel(1,4),0.0_dp,1.0e-12_dp,'group centered 2')
  print '(a)','test_cross_sectional: PASS'
contains
  subroutine assert_true(ok,msg)
    logical,intent(in)::ok;character(len=*),intent(in)::msg
    if(.not.ok)then;write(*,'(a)')'FAIL: '//msg;error stop 1;end if
  end subroutine
  subroutine assert_close(a,b,tol,msg)
    real(dp),intent(in)::a,b,tol;character(len=*),intent(in)::msg
    call assert_true(abs(a-b)<=tol,msg)
  end subroutine
end program test_cross_sectional
