program demo_portvine
   use portvine
   use rugarch, only : garch_spec,dist_norm,seed_rng
   use rvinecopulib, only : bicop_indep,bicop_gaussian
   implicit none
   integer,parameter::d=3,n=70
   real(dp)::returns(d,n),weights(d,4),u1,u2,common
   type(garch_spec)::spec
   type(marginal_settings_type)::marginal_settings
   type(vine_settings_type)::vine_settings
   type(portvine_roll_result)::result
   integer::i,status,last

   call seed_rng(20260802)
   do i=1,n
      call random_number(u1);call random_number(u2)
      common=0.015_dp*(u1-0.5_dp)
      returns(1,i)=common+0.004_dp*(u2-0.5_dp)
      call random_number(u1);returns(2,i)=0.7_dp*common+0.008_dp*(u1-0.5_dp)
      call random_number(u1);returns(3,i)=0.010_dp*(u1-0.5_dp)
   end do
   spec=make_portvine_spec(0,0,1,1,cond_dist=dist_norm)
   marginal_settings=make_marginal_settings(50,10,d,spec)
   marginal_settings%max_iterations=140
   vine_settings=make_vine_settings(30,5,vine_dvine,[bicop_indep,bicop_gaussian])
   weights=0.0_dp
   weights(2,:)=0.60_dp;weights(3,:)=0.40_dp
   call estimate_risk_roll(returns,marginal_settings,vine_settings,[0.01_dp,0.05_dp], &
      [risk_var,risk_es_mean],80,result,weights=weights,cond_indices=[1], &
      cond_u=[0.01_dp,0.05_dp,0.50_dp],n_mc_samples=100, &
      prior_residual_strategy=.true.,status=status)
   if(status/=0)then
      print '(a)',trim(result%message)
      stop 1
   end if
   last=size(result%realized)
   print '(a,i0)','forecast observations: ',last
   print '(a,2f11.6)','last overall VaR:      ',result%overall(1,:,last)
   print '(a,2f11.6)','last stress VaR (1%):  ',result%conditional(1,:,1,last)
   print '(a,f11.6)','last realized return:  ',result%realized(last)
end program demo_portvine
