program ffp_demo
  use ffp_mod
  implicit none
  integer,parameter::n=7
  real(dp)::x(n,1),prior(n),aeq(2,n),beq(2),post(n),mu(1),sigma(1,1)
  integer::info,i
  x(:,1)=[-3.0_dp,-2.0_dp,-1.0_dp,0.0_dp,1.0_dp,2.0_dp,3.0_dp]
  prior=1.0_dp/real(n,dp); aeq(1,:)=1.0_dp; beq(1)=1.0_dp; aeq(2,:)=x(:,1); beq(2)=0.75_dp
  call entropy_pool_equalities(prior,aeq,beq,post,info)
  call weighted_moments(x,post,mu,sigma)
  print '(a,i0)', 'solver info: ',info
  print '(a,f10.6)', 'posterior mean: ',mu(1)
  print '(a,f10.6)', 'effective scenarios: ',effective_scenarios(post)
  print '(a)', 'scenario  probability'
  do i=1,n; print '(f8.3,2x,f12.8)',x(i,1),post(i); end do
end program ffp_demo
