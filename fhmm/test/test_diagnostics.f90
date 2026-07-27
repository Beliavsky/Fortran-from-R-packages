! SPDX-License-Identifier: GPL-3.0-only
program test_diagnostics
   use fhmm
   implicit none
   type(hmm_parameters) :: par,reordered
   type(forecast_result) :: fc,fc_upstream
   type(model_comparison) :: mc
   real(dp) :: probs(2)
   integer, allocatable :: chunks(:)

   par%distribution=dist_normal
   allocate(par%gamma(2,2),par%mu(2),par%sigma(2),par%df(2))
   par%gamma=reshape([0.9_dp,0.2_dp,0.1_dp,0.8_dp],[2,2])
   par%mu=[1.0_dp,-1.0_dp];par%sigma=[0.5_dp,0.5_dp];par%df=10.0_dp
   reordered=reorder_hmm_states(par)
   if(maxval(abs(reordered%mu-[-1.0_dp,1.0_dp]))>1.0e-15_dp)error stop 1
   probs=[1.0_dp,0.0_dp]
   fc=forecast_hmm(par,3,alpha=0.05_dp,last_probabilities=probs)
   fc_upstream=forecast_hmm(par,3,alpha=0.05_dp,last_probabilities=probs,upstream_quantiles=.true.)
   if(.not.fc%ok.or..not.fc_upstream%ok)error stop 1
   if(maxval(abs(sum(fc%state_probabilities,dim=1)-1.0_dp))>1.0e-14_dp)error stop 1
   if(any(fc%lower>fc%median).or.any(fc%median>fc%upper))error stop 1
   mc=compare_hmm_model(-100.0_dp,6,1000)
   call assert_close(mc%aic,212.0_dp,1.0e-14_dp,'AIC')
   call assert_close(mc%bic,200.0_dp+6.0_dp*log(1000.0_dp),1.0e-14_dp,'BIC')
   chunks=compute_chunk_lengths(20,-1,'w',seed=1)
   if(any(chunks<1).or.any(chunks>5))error stop 1
   if(any(compute_chunk_lengths(4,7,'m')/=7))error stop 1
   print '(a)','test_diagnostics: PASS'
contains
   subroutine assert_close(x,y,tol,label)
      real(dp),intent(in)::x,y,tol;character(len=*),intent(in)::label
      if(abs(x-y)>tol)then;write(*,*)trim(label),x,y;error stop 1;end if
   end subroutine
end program test_diagnostics
