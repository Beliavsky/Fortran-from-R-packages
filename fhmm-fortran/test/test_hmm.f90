! SPDX-License-Identifier: GPL-3.0-only
program test_hmm
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use fhmm
   implicit none
   type(hmm_parameters) :: par
   type(inference_result) :: inf
   type(hmm_simulation) :: sim1,sim2
   integer, allocatable :: states(:)
   real(dp), allocatable :: residuals(:)
   real(dp) :: obs(4),ll

   par%distribution=dist_normal
   allocate(par%gamma(2,2),par%mu(2),par%sigma(2),par%df(2))
   par%gamma=reshape([0.9_dp,0.2_dp,0.1_dp,0.8_dp],[2,2])
   par%mu=[-1.0_dp,1.0_dp];par%sigma=0.5_dp;par%df=10.0_dp
   obs=[-1.0_dp,-0.5_dp,1.0_dp,1.5_dp]
   ll=hmm_log_likelihood(obs,par)
   call assert_close(ll,-4.923139050080705_dp,2.0e-13_dp,'log likelihood')
   inf=forward_backward(obs,par)
   if(.not.inf%ok)then;write(*,*)trim(inf%message);error stop 1;end if
   call assert_close(inf%log_likelihood,ll,2.0e-13_dp,'forward likelihood')
   if(maxval(abs(sum(inf%filtered,dim=1)-1.0_dp))>2.0e-14_dp)error stop 1
   if(maxval(abs(sum(inf%smoothed,dim=1)-1.0_dp))>2.0e-14_dp)error stop 1
   states=viterbi_decode(obs,par)
   if(any(states/=[1,1,2,2]))then;write(*,*)states;error stop 1;end if
   residuals=pseudo_residuals(obs,states,par)
   if(any(.not.ieee_is_finite(residuals)))error stop 1
   sim1=simulate_hmm_model(par,100,seed=731)
   sim2=simulate_hmm_model(par,100,seed=731)
   if(any(sim1%states/=sim2%states))error stop 1
   if(maxval(abs(sim1%observations-sim2%observations))>0.0_dp)error stop 1
   print '(a)','test_hmm: PASS'
contains
   subroutine assert_close(x,y,tol,label)
      real(dp),intent(in)::x,y,tol;character(len=*),intent(in)::label
      if(abs(x-y)>tol)then;write(*,*)trim(label),x,y,abs(x-y);error stop 1;end if
   end subroutine
end program test_hmm
