! SPDX-License-Identifier: GPL-3.0-only
program test_views
  use portfolio_analytics
  use test_support
  implicit none
  real(dp) :: mu(2),sigma(2,2),pick(1,2),views(1),pmu(2),psig(2,2)
  real(dp) :: prior(4),aeq(2,4),beq(2),returns(4,2),wmu(2),wsig(2,2),posterior(4)
  real(dp) :: ranking_mu(2),ranking_sigma(2,2),cent(4),acmu(2)
  integer :: info,order2(2)
  logical :: converged
  type(entropy_result) :: ep

  mu=[0.02_dp,0.01_dp]
  sigma=reshape([0.04_dp,0.01_dp,0.01_dp,0.09_dp],[2,2])
  pick(1,:)=[1.0_dp,-1.0_dp]
  views=[0.05_dp]
  call black_litterman(mu,sigma,pick,views,pmu,psig,info=info)
  call assert_true(info==0,'Black-Litterman solve')
  call assert_true(dot_product(pick(1,:),pmu)>dot_product(pick(1,:),mu),'view shifts posterior')
  call assert_true(all([(psig(info,info)>0.0_dp,info=1,2)]),'posterior covariance diagonal')

  prior=0.25_dp
  aeq(1,:)=1.0_dp
  aeq(2,:)=[0.0_dp,1.0_dp,2.0_dp,3.0_dp]
  beq=[1.0_dp,2.0_dp]
  call entropy_pool(prior,ep,aeq=aeq,beq=beq)
  call assert_true(ep%converged,'entropy equality convergence')
  call assert_close(sum(ep%probabilities),1.0_dp,1.0e-10_dp,'posterior probability sum')
  call assert_close(dot_product(ep%probabilities,aeq(2,:)),2.0_dp,1.0e-10_dp,'posterior mean view')

  returns(:,1)=[-0.02_dp,-0.01_dp,0.01_dp,0.02_dp]
  returns(:,2)=[0.02_dp,0.01_dp,-0.01_dp,-0.02_dp]
  call meucci_moments(returns,prior,wmu,wsig)
  call assert_all_close(wmu,[0.0_dp,0.0_dp],1.0e-14_dp,'Meucci prior moments')
  order2=[1,2]
  call meucci_ranking(returns,prior,order2,ranking_mu,ranking_sigma,posterior,converged)
  call assert_true(converged,'Meucci ranking convergence')
  call assert_true(ranking_mu(1)<=ranking_mu(2)+1.0e-9_dp,'ranking expected return order')
  call centroid(4,cent)
  call assert_true(all(cent(1:3)>cent(2:4)),'centroid descending')
  call ac_ranking(returns,order2,acmu)
  call assert_true(acmu(1)<acmu(2),'AC ranking order')
  print '(a)','test_views: PASS'
end program test_views
