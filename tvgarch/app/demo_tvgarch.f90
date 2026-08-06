program demo_tvgarch
   use tvgarch, only : dp, tvgarch_spec, tvgarch_simulation, tvgarch_fit, &
                       make_tvgarch_spec, make_tv_spec, tvgarch_simulate, fit_tvgarch, &
                       tvgarch_forecast, set_random_seed
   implicit none
   integer,parameter::n=400
   type(tvgarch_spec)::spec
   type(tvgarch_simulation)::sim
   type(tvgarch_fit)::fit
   real(dp)::ph(3)
   real(dp),allocatable::x(:),forecast(:)
   integer::i,st
   allocate(x(n));do i=1,n;x(i)=real(i,dp)/real(n,dp);end do
   call make_tvgarch_spec(spec,order_g=[1],order_h=[1,1,0],speed_option=2,status=st)
   call make_tv_spec(spec%tv,[1],1.0_dp,[2.0_dp],[log(20.0_dp)],[0.55_dp],2)
   ph=[0.10_dp,0.08_dp,0.85_dp]
   call set_random_seed(20260804)
   call tvgarch_simulate(n,spec,ph,sim,xtv=x)
   call fit_tvgarch(sim%y,spec,fit,xtv=x,initial_g=[1.0_dp,1.5_dp,log(15.0_dp),0.5_dp], &
                    initial_h=ph,max_outer=10,turbo=.true.)
   call tvgarch_forecast(fit,5,forecast,st,n_sim=250)
   write(*,'(a,4f10.4)')'TV estimates: ',fit%par_g
   write(*,'(a,3f10.4)')'GARCH estimates: ',fit%hfit%par
   write(*,'(a,f12.4)')'Log likelihood: ',fit%loglik
   write(*,'(a,5f10.4)')'Variance forecast: ',forecast
end program demo_tvgarch
