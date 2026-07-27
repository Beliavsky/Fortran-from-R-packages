! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of MSGARCH, copyright (C) MSGARCH authors.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
program test_filter_simulation_risk
   use msgarch
   use test_helpers
   implicit none
   type(msgarch_spec)::spec
   type(simulation_result)::simulation,ahead
   type(filter_result)::filtered
   type(risk_result)::risk,itsrisk
   real(dp),allocatable::pit(:),vol(:),pdf(:,:),cdf(:,:),state_forecast(:,:),power(:,:),vf(:),mf(:)
   real(dp),allocatable::forecast_pdf(:,:),forecast_cdf(:,:),forecast_draw(:,:)
   real(dp)::p1,p2,uncvol
   integer::t

   call seed_rng(314159)
   spec=create_spec([character(len=12)::'sGARCH','gjrGARCH'],[character(len=8)::'norm','std'])
   spec%regime(1)%omega=0.04_dp;spec%regime(1)%alpha=0.08_dp;spec%regime(1)%beta=0.88_dp
   spec%regime(2)%omega=0.12_dp;spec%regime(2)%alpha=0.10_dp;spec%regime(2)%gamma=0.08_dp
   spec%regime(2)%beta=0.72_dp;spec%regime(2)%shape=9.0_dp
   spec%transition=reshape([0.96_dp,0.04_dp,0.08_dp,0.92_dp],[2,2],order=[2,1])
   call assert_true(spec_valid(spec),'two-regime specification valid')
   simulation=simulate_msgarch(spec,500,3)
   filtered=hamilton_filter(spec,simulation%draw(1,:))
   call assert_all_finite(filtered%log_density,'state log densities finite')
   call assert_true(filtered%loglik>-1.0e100_dp,'Hamilton log likelihood finite')
   do t=1,500
      call assert_close(sum(filtered%predicted(t,:)),1.0_dp,1.0e-10_dp,'predicted probabilities sum one')
      call assert_close(sum(filtered%filtered(t,:)),1.0_dp,1.0e-10_dp,'filtered probabilities sum one')
      call assert_close(sum(filtered%smoothed(t,:)),1.0_dp,1.0e-9_dp,'smoothed probabilities sum one')
   end do
   call assert_true(all(filtered%viterbi>=1).and.all(filtered%viterbi<=2),'Viterbi states in range')
   call assert_close(sum(filtered%next_probability),1.0_dp,1.0e-12_dp,'next probabilities sum one')

   p1=predictive_cdf(spec,simulation%draw(1,:),-1.0_dp)
   p2=predictive_cdf(spec,simulation%draw(1,:),1.0_dp)
   call assert_true(p1<p2,'predictive CDF monotone')
   call assert_true(predictive_pdf(spec,simulation%draw(1,:),0.0_dp)>0.0_dp,'predictive density positive')
   pit=pit_values(spec,simulation%draw(1,:));call assert_true(all(pit>=0.0_dp).and.all(pit<=1.0_dp),'PIT in unit interval')
   vol=conditional_volatility(spec,simulation%draw(1,:));call assert_true(all(vol>0.0_dp),'conditional volatility positive')
   pdf=in_sample_pdf(spec,simulation%draw(1,:),[-1.0_dp,0.0_dp,1.0_dp])
   cdf=in_sample_cdf(spec,simulation%draw(1,:),[-1.0_dp,0.0_dp,1.0_dp])
   call assert_true(all(pdf>=0.0_dp),'in-sample densities nonnegative')
   call assert_true(all(cdf>=0.0_dp).and.all(cdf<=1.0_dp),'in-sample CDF in range')

   ahead=simulate_ahead(spec,simulation%draw(1,:),6,1200)
   call assert_all_finite(ahead%draw,'ahead simulations finite')
   vf=forecast_volatility(spec,simulation%draw(1,:),6,1200)
   mf=forecast_mean(spec,simulation%draw(1,:),6,1200)
   call assert_true(all(vf>0.0_dp),'volatility forecasts positive')
   call assert_true(maxval(abs(mf))<0.5_dp,'forecast means near zero')
   call predictive_distribution_forecast(spec,simulation%draw(1,:),[-1.0_dp,0.0_dp,1.0_dp],4,1600, &
      forecast_pdf,forecast_cdf,forecast_draw)
   call assert_true(all(forecast_pdf>=0.0_dp),'predictive density forecasts nonnegative')
   call assert_true(all(forecast_cdf>=0.0_dp).and.all(forecast_cdf<=1.0_dp),'predictive CDF forecasts in range')
   call assert_true(all(forecast_cdf(:,1)<=forecast_cdf(:,2)).and.all(forecast_cdf(:,2)<=forecast_cdf(:,3)), &
      'predictive CDF forecasts monotone')
   call assert_close(forecast_pdf(1,2),predictive_pdf(spec,simulation%draw(1,:),0.0_dp),1.0e-12_dp, &
      'one-step analytical predictive density')
   call assert_true(all(shape(forecast_draw)==[1600,4]),'predictive simulations retained')
   uncvol=unconditional_volatility(spec,nsim=120,n_ahead=120,burn=40)
   call assert_true(uncvol>0.0_dp.and.uncvol<20.0_dp,'simulation-based unconditional volatility')

   state_forecast=forecast_state_probabilities(spec,simulation%draw(1,:),8)
   do t=1,8;call assert_close(sum(state_forecast(t,:)),1.0_dp,1.0e-12_dp,'forecast state probabilities');end do
   power=transition_matrix_power(spec%transition,5)
   call assert_true(all(abs(sum(power,dim=2)-1.0_dp)<1.0e-12_dp),'transition matrix power stochastic')

   risk=risk_forecast(spec,simulation%draw(1,:),[0.01_dp,0.05_dp],4,5000,return_simulations=.true.)
   call assert_true(all(risk%es<=risk%var+1.0e-12_dp),'expected shortfall below VaR')
   call assert_true(size(risk%simulations,1)==5000,'risk simulations retained')
   itsrisk=risk_in_sample(spec,simulation%draw(1,1:80),[0.05_dp])
   call assert_true(all(itsrisk%es(2:,1)<=itsrisk%var(2:,1)+1.0e-8_dp),'in-sample ES below VaR')

   write(*,'(a)')'Filtering, simulation, forecast, PIT, and risk tests passed.'
end program test_filter_simulation_risk
