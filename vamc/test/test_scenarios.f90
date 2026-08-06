program test_scenarios
  use vamc
  use vamc_test_support
  implicit none
  real(dp) :: covariance(2,2), forward(6), fund_map(3,2)
  real(dp), allocatable :: x1(:,:,:), x2(:,:,:), funds(:,:,:), two_d(:,:)
  type(status_type) :: status
  covariance=reshape([0.04_dp,0.01_dp,0.01_dp,0.09_dp],[2,2])
  forward=0.02_dp
  call gen_index_scen(covariance,4,6,1.0_dp/12.0_dp,forward,x1,seed=123,status=status)
  call assert_true(status%ok(),'generate index scenarios')
  call gen_index_scen(covariance,4,6,1.0_dp/12.0_dp,forward,x2,seed=123,status=status)
  call assert_true(maxval(abs(x1-x2))<tiny(1.0_dp),'seeded scenario reproducibility')
  call assert_true(all(x1>0.0_dp),'positive Black-Scholes factors')
  fund_map(1,:)=[1.0_dp,0.0_dp]
  fund_map(2,:)=[0.0_dp,1.0_dp]
  fund_map(3,:)=[0.6_dp,0.4_dp]
  call gen_fund_scen(fund_map,x1,funds,status)
  call assert_true(status%ok(),'generate fund scenarios')
  call assert_close(funds(1,1,3),0.6_dp*x1(1,1,1)+0.4_dp*x1(1,1,2),1.0e-14_dp,'fund mapping')
  call gen_fund_scen(fund_map,x1(1,:,:),two_d,status)
  call assert_true(status%ok(),'two-dimensional fund mapping')
  call assert_all_close(two_d(:,3),funds(1,:,3),1.0e-14_dp,'rank-specific fund mapping')
  print '(a)','test_scenarios: PASS'
end program test_scenarios
