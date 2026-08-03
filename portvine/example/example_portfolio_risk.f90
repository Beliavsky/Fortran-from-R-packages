program example_portfolio_risk
   use portvine
   use rugarch, only : garch_spec,dist_norm,seed_rng
   use rvinecopulib, only : bicop_indep,bicop_gaussian
   implicit none
   integer,parameter::d=3,n=70
   real(dp)::r(d,n),w(d,4),u1,u2
   type(garch_spec)::spec
   type(marginal_settings_type)::ms
   type(vine_settings_type)::vs
   type(portvine_roll_result)::ans
   integer::i,status
   call seed_rng(271828)
   do i=1,n
      call random_number(u1);call random_number(u2)
      r(1,i)=0.012_dp*(u1-0.5_dp)
      r(2,i)=0.7_dp*r(1,i)+0.007_dp*(u2-0.5_dp)
      call random_number(u1);r(3,i)=0.010_dp*(u1-0.5_dp)
   end do
   spec=make_portvine_spec(0,0,1,1,cond_dist=dist_norm)
   ms=make_marginal_settings(50,10,d,spec);ms%max_iterations=140
   vs=make_vine_settings(30,5,vine_dvine,[bicop_indep,bicop_gaussian])
   w=0.0_dp;w(2,:)=0.6_dp;w(3,:)=0.4_dp
   call estimate_risk_roll(r,ms,vs,[0.05_dp],[risk_var,risk_es_mean],100,ans, &
      weights=w,cond_indices=[1],cond_u=[0.05_dp],prior_residual_strategy=.true.,status=status)
   print '(a,f10.6)','last conditional VaR: ',ans%conditional(1,1,1,size(ans%realized))
   print '(a,f10.6)','last conditional ES:  ',ans%conditional(2,1,1,size(ans%realized))
end program example_portfolio_risk
