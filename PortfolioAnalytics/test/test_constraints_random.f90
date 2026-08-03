! SPDX-License-Identifier: GPL-3.0-only
program test_constraints_random
  use portfolio_analytics
  use test_support
  implicit none
  type(portfolio_constraints) :: c
  real(dp) :: w(4),candidate(4),groups(2),mu(4)
  real(dp),allocatable :: samples(:,:),grid(:,:)
  logical :: ok
  integer :: i,nfound

  call initialize_constraints(c,4,max_weight=[0.6_dp,0.6_dp,0.6_dp,0.6_dp])
  allocate(c%group_a(2,4),c%group_lower(2),c%group_upper(2))
  c%group_a=0.0_dp
  c%group_a(1,1:2)=1.0_dp
  c%group_a(2,3:4)=1.0_dp
  c%group_lower=[0.2_dp,0.2_dp]
  c%group_upper=[0.8_dp,0.8_dp]
  candidate=[2.0_dp,-1.0_dp,0.1_dp,0.2_dp]
  call repair_weights(candidate,c,w,ok)
  call assert_true(ok,'box-sum repair')
  call assert_close(sum(w),1.0_dp,1.0e-10_dp,'repaired sum')
  mu=[0.01_dp,0.02_dp,0.03_dp,0.04_dp]
  allocate(samples(4,40),grid(4,20))
  call random_portfolios(c,40,samples,nfound,777,mu)
  call assert_true(nfound==40,'random portfolio count')
  do i=1,nfound
    call assert_true(check_constraints(samples(:,i),c,mu),'random portfolio feasibility')
    call group_exposures(samples(:,i),c,groups)
    call assert_true(all(groups>=0.2_dp-1.0e-8_dp),'group lower')
    call assert_true(all(groups<=0.8_dp+1.0e-8_dp),'group upper')
  end do
  call random_grid_portfolios(c,20,0.05_dp,grid,nfound,778,mu)
  call assert_true(nfound>10,'grid portfolios generated')
  c%max_positions=2
  candidate=[0.4_dp,0.3_dp,0.2_dp,0.1_dp]
  call repair_weights(candidate,c,w,ok)
  call assert_true(ok,'position repair')
  call assert_true(count_positions(w)<=2,'position limit')
  print '(a)','test_constraints_random: PASS'
end program test_constraints_random
