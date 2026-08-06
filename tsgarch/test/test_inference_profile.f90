program test_inference_profile
  use ghyp_kinds, only : dp, i8
  use tsgarch
  use test_support
  implicit none
  real(dp)::dummy(30),grid(5)
  real(dp),allocatable::pit(:),lo(:),hi(:),cov(:,:)
  type(garch_spec)::spec
  type(garch_parameters)::par
  type(garch_simulation)::sim
  type(garch_fit)::fit
  type(profile_result)::prof
  integer::i,status
  dummy=[(0.01_dp*sin(real(i,dp)),i=1,30)]
  spec=standard_spec('garch','norm')
  par=standard_parameters(dummy,spec)
  sim=simulate_garch(spec,par,100,burn=20,seed=19_i8)
  fit%spec=spec
  fit%parameters=par
  fit%filtered=filter_garch(sim%series(:,1),spec,par)
  fit%status=tsg_success
  call pack_parameters(spec,par,fit%packed_parameters,fit%parameter_names)
  fit%npars=size(fit%packed_parameters)
  fit%log_likelihood=fit%filtered%log_likelihood
  call numerical_inference(sim%series(:,1),spec,par,fit)
  call assert_true(allocated(fit%hessian).and.allocated(fit%scores),'inference arrays absent')
  pit=probability_integral_transform(fit)
  call assert_true(all(pit>=0.0_dp).and.all(pit<=1.0_dp),'PIT outside unit interval')
  call confidence_intervals(fit,0.95_dp,lo,hi,status)
  call assert_true(status==tsg_success.and.all(hi>=lo),'confidence intervals failed')
  call covariance_opg(fit,cov,status)
  call assert_true(status==tsg_success.or.status==tsg_singular,'OPG status')
  grid=par%omega*[0.7_dp,0.85_dp,1.0_dp,1.15_dp,1.3_dp]
  prof=profile_likelihood(sim%series(:,1),fit,2,grid)
  call assert_true(prof%status==tsg_success.and.all(prof%log_likelihood>-huge(1.0_dp)),'profile failed')
  call assert_true(half_life(fit)>0.0_dp,'half life failed')
  write(*,'(a)')'test_inference_profile: PASS'
end program test_inference_profile
