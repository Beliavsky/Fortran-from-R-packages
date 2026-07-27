! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of MSGARCH, copyright (C) MSGARCH authors.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
program test_estimation_hmm
   use msgarch
   use test_helpers
   implicit none
   type(msgarch_spec)::true_spec,mix_spec,single_spec,mean_spec
   type(simulation_result)::simulation,mix_sim
   type(filter_result)::start_filter
   type(fit_result)::fit,fixed_fit,tied_fit
   type(mcmc_result)::chain
   type(posterior_state_result)::posterior_state
   type(risk_result)::posterior_risk
   type(hmm_fit_result)::hmm,mix
   real(dp),allocatable::start(:),proposal(:),hmm_y(:),post_vol(:),post_pit(:)
   logical,allocatable::fixed(:)
   integer,allocatable::tie(:),states(:)
   integer::i,n
   real(dp)::post_uncvol

   call seed_rng(86420)
   true_spec=create_spec([character(len=12)::'sGARCH'],[character(len=8)::'norm'])
   true_spec%regime(1)%omega=0.06_dp;true_spec%regime(1)%alpha=0.10_dp;true_spec%regime(1)%beta=0.84_dp
   simulation=simulate_msgarch(true_spec,420,1)
   start=pack_parameters(true_spec);start=[0.10_dp,0.16_dp,0.70_dp]
   call unpack_parameters(true_spec,start,mix_spec)
   start_filter=hamilton_filter(mix_spec,simulation%draw(1,:))
   fit=fit_ml(true_spec,simulation%draw(1,:),start=start,max_iterations=500,tolerance=2.0e-7_dp)
   call assert_true(fit%loglik>=start_filter%loglik-1.0e-6_dp,'ML improves likelihood')
   call assert_true(all(fit%parameters>0.0_dp),'ML parameters positive')
   call assert_all_finite(fit%hessian,'ML Hessian finite')
   call assert_all_finite(fit%covariance,'ML covariance finite')
   call assert_true(fit%aic<huge(1.0_dp)/2.0_dp.and.fit%bic<huge(1.0_dp)/2.0_dp,'AIC and BIC finite')

   allocate(fixed(3));fixed=.true.
   fixed_fit=fit_ml(true_spec,simulation%draw(1,:),start=pack_parameters(true_spec),fixed_mask=fixed, &
      fixed_values=pack_parameters(true_spec))
   call assert_true(fixed_fit%converged,'all-fixed evaluation succeeds')
   call assert_close(maxval(abs(fixed_fit%parameters-pack_parameters(true_spec))),0.0_dp,1.0e-14_dp,'fixed parameters preserved')

   mix_spec=create_spec([character(len=12)::'sARCH','sARCH'],[character(len=8)::'norm','norm'],.true.)
   mix_spec%regime(1)%omega=0.05_dp;mix_spec%regime(1)%alpha=0.12_dp
   mix_spec%regime(2)%omega=0.20_dp;mix_spec%regime(2)%alpha=0.30_dp
   mix_spec%transition=0.0_dp;mix_spec%transition(:,1)=0.65_dp;mix_spec%transition(:,2)=0.35_dp
   mix_sim=simulate_msgarch(mix_spec,350,1)
   start=pack_parameters(mix_spec);allocate(tie(size(start)));tie=0
   tie(1)=1;tie(3)=1
   tied_fit=fit_ml(mix_spec,mix_sim%draw(1,:),start=start,tie_group=tie,max_iterations=300)
   single_spec=extract_regime(mix_spec,2)
   call assert_true(single_spec%k==1,'single-regime extraction count')
   call assert_true(spec_valid(single_spec),'single-regime extraction validity')
   call assert_close(tied_fit%parameters(1),tied_fit%parameters(3),1.0e-12_dp,'regime-constant tied parameter')
   call assert_true(tied_fit%loglik>-1.0e100_dp,'tied mixture fit finite')

   proposal=[0.008_dp,0.008_dp,0.008_dp]
   chain=fit_mcmc(true_spec,simulation%draw(1,1:260),700,200,5,start=pack_parameters(true_spec),proposal_sd=proposal)
   call assert_true(size(chain%draws,1)==100,'MCMC burn-in and thinning')
   call assert_true(chain%acceptance_rate>0.01_dp.and.chain%acceptance_rate<0.99_dp,'MCMC chain moves')
   call assert_all_finite(chain%draws,'MCMC draws finite')
   call assert_all_finite(chain%posterior_mean,'posterior means finite')
   call assert_true(chain%dic<huge(1.0_dp)/2.0_dp,'DIC finite')
   mean_spec=posterior_mean_spec(true_spec,chain%draws)
   call assert_true(spec_valid(mean_spec),'posterior mean specification valid')
   posterior_state=posterior_state_probabilities(true_spec,simulation%draw(1,1:120),chain%draws(1:20,:))
   call assert_true(all(abs(sum(posterior_state%filtered,dim=2)-1.0_dp)<1.0e-10_dp),'posterior filtered probabilities')
   call assert_true(posterior_predictive_pdf(true_spec,simulation%draw(1,1:120),chain%draws(1:20,:),0.0_dp)>0.0_dp, &
      'posterior predictive density')
   call assert_true(posterior_predictive_cdf(true_spec,simulation%draw(1,1:120),chain%draws(1:20,:),0.0_dp)>0.0_dp, &
      'posterior predictive CDF')
   post_vol=posterior_volatility(true_spec,simulation%draw(1,1:120),chain%draws(1:20,:))
   post_pit=posterior_pit(true_spec,simulation%draw(1,1:120),chain%draws(1:20,:))
   call assert_true(all(post_vol>0.0_dp),'posterior volatility positive')
   call assert_true(all(post_pit>=0.0_dp).and.all(post_pit<=1.0_dp),'posterior PIT in range')
   post_uncvol=posterior_unconditional_volatility(true_spec,chain%draws(1:5,:),nsim=60,n_ahead=80,burn=30)
   call assert_true(post_uncvol>0.0_dp.and.post_uncvol<20.0_dp,'posterior unconditional volatility')
   posterior_risk=posterior_risk_forecast(true_spec,chain%draws(1:20,:),simulation%draw(1,1:120),[0.05_dp],2,20)
   call assert_true(all(posterior_risk%es<=posterior_risk%var),'posterior expected shortfall')

   n=500;allocate(hmm_y(n))
   do i=1,250;hmm_y(i)=-0.8_dp+0.45_dp*random_normal();end do
   do i=251,n;hmm_y(i)=0.9_dp+0.65_dp*random_normal();end do
   hmm=fit_gaussian_hmm(hmm_y,2,max_iterations=300)
   call assert_true(hmm%loglik>-1.0e100_dp,'Gaussian HMM EM finite')
   call assert_true(all(hmm%variance>0.0_dp),'HMM variances positive')
   call assert_true(all(abs(sum(hmm%transition,dim=2)-1.0_dp)<1.0e-10_dp),'HMM transition rows sum one')
   states=viterbi_gaussian_hmm(hmm_y,hmm%mean,hmm%variance,hmm%transition)
   call assert_true(all(states>=1).and.all(states<=2),'Gaussian HMM Viterbi states')
   mix=fit_gaussian_mixture(hmm_y,2,max_iterations=300)
   call assert_close(sum(mix%probability),1.0_dp,1.0e-10_dp,'Gaussian mixture weights sum one')
   call assert_true(all(mix%variance>0.0_dp),'Gaussian mixture variances positive')

   write(*,'(a)')'ML, constrained fitting, MCMC, DIC, and Gaussian EM tests passed.'
end program test_estimation_hmm
