program mean_view
  use ffp_mod
  implicit none
  real(dp)::x(5,1),p0(5),p(5),a(2,5),b(2),mu(1),s(1,1)
  integer::info
  x(:,1)=[-2.0_dp,-1.0_dp,0.0_dp,1.0_dp,2.0_dp]; p0=0.2_dp
  a(1,:)=1.0_dp; b(1)=1.0_dp; a(2,:)=x(:,1); b(2)=0.5_dp
  call entropy_pool_equalities(p0,a,b,p,info)
  call weighted_moments(x,p,mu,s)
  print '(a,f10.6)', 'target/realized mean: ',mu(1)
end program mean_view
