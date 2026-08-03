! SPDX-License-Identifier: MIT
! PortfolioTesteR modern Fortran translation
program example_walk_forward
  use portfolio_tester
  implicit none
  real(dp),allocatable::prices(:,:)
  real(dp)::params(2,6)
  type(walk_forward_result)::wf
  call generate_sample_prices(220,10,prices,4401_i8)
  params(:,1)=[8.0_dp,3.0_dp];params(:,2)=[12.0_dp,3.0_dp]
  params(:,3)=[20.0_dp,3.0_dp];params(:,4)=[8.0_dp,5.0_dp]
  params(:,5)=[12.0_dp,5.0_dp];params(:,6)=[20.0_dp,5.0_dp]
  call run_walk_forward(prices,params,momentum_top_n_strategy,104,26,26,wf, &
    warmup_periods=20,cost_bps=5.0_dp,frequency=52.0_dp)
  print '(a,i0)','walk-forward windows: ',size(wf%best_index)
  print '(a,6i4)','chosen grid indices: ',wf%best_index
  if(size(wf%stitched_value)>0)then
    print '(a,f10.4)','stitched OOS return: ', &
      wf%stitched_value(size(wf%stitched_value))/wf%stitched_value(1)-1.0_dp
  end if
end program example_walk_forward
