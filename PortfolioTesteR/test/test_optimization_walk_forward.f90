! SPDX-License-Identifier: MIT
! PortfolioTesteR modern Fortran translation
program test_optimization_walk_forward
  use portfolio_tester
  implicit none
  real(dp),allocatable::prices(:,:)
  real(dp)::params(2,4)
  type(grid_result)::grid
  type(walk_forward_result)::wf
  call generate_sample_prices(140,6,prices,7123_i8)
  params(:,1)=[4.0_dp,2.0_dp];params(:,2)=[8.0_dp,2.0_dp]
  params(:,3)=[4.0_dp,3.0_dp];params(:,4)=[8.0_dp,3.0_dp]
  call run_param_grid(prices,params,momentum_top_n_strategy,grid)
  call assert_true(grid%best_index>=1.and.grid%best_index<=4,'grid best index')
  call run_walk_forward(prices,params,momentum_top_n_strategy,60,20,20,wf,warmup_periods=10)
  call assert_true(size(wf%best_index)>=2,'walk forward windows')
  call assert_true(count(wf%best_index>0)>0,'walk forward selections')
  call assert_true(size(wf%stitched_value)>1,'stitched equity')
  print '(a)','test_optimization_walk_forward: PASS'
contains
  subroutine assert_true(ok,msg)
    logical,intent(in)::ok;character(len=*),intent(in)::msg
    if(.not.ok)then;write(*,'(a)')'FAIL: '//msg;error stop 1;end if
  end subroutine
end program test_optimization_walk_forward
