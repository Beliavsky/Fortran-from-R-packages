! SPDX-License-Identifier: GPL-3.0-only
program test_hhmm
   use fhmm
   implicit none
   type(hhmm_parameters) :: par
   type(hhmm_inference_result) :: inf
   type(hhmm_simulation) :: sim
   type(hhmm_fit_result) :: fit
   type(fit_options) :: opt
   integer, allocatable :: cs(:),fs(:,:)
   real(dp) :: coarse(2),fine(2,3),ll
   integer :: lengths(2)

   call make_model(par)
   coarse=[-1.0_dp,1.0_dp]
   fine=0.0_dp;fine(1,1:2)=[-2.0_dp,-1.8_dp];fine(2,1:2)=[2.0_dp,1.8_dp]
   lengths=[2,2]
   ll=hhmm_log_likelihood(coarse,fine,lengths,par)
   call assert_close(ll,-5.221151151411109_dp,3.0e-12_dp,'HHMM likelihood')
   inf=hhmm_forward_backward(coarse,fine,lengths,par)
   if(.not.inf%ok)then;write(*,*)trim(inf%message);error stop 1;end if
   if(maxval(abs(sum(inf%coarse_smoothed,dim=1)-1.0_dp))>1.0e-12_dp)error stop 1
   call decode_hhmm(coarse,fine,lengths,par,cs,fs)
   if(any(cs/=[1,2]))then;write(*,*)cs;error stop 1;end if
   sim=simulate_hhmm_model(par,[3,2,4],seed=8)
   if(any(sim%chunk_lengths/=[3,2,4]))error stop 1
   if(size(sim%fine_observations,2)/=4)error stop 1
   opt%runs=1;opt%max_iterations=5;opt%compute_hessian=.false.
   fit=fit_hhmm(coarse,fine,lengths,2,2,dist_normal,dist_normal,opt,par)
   if(.not.fit%ok)then;write(*,*)trim(fit%message);error stop 1;end if
   print '(a)','test_hhmm: PASS'
contains
   subroutine make_model(p)
      type(hhmm_parameters),intent(out)::p
      integer::s
      p%coarse%distribution=dist_normal
      allocate(p%coarse%gamma(2,2),p%coarse%mu(2),p%coarse%sigma(2),p%coarse%df(2),p%fine(2))
      p%coarse%gamma=reshape([0.9_dp,0.2_dp,0.1_dp,0.8_dp],[2,2])
      p%coarse%mu=[-1.0_dp,1.0_dp];p%coarse%sigma=0.5_dp;p%coarse%df=10.0_dp
      do s=1,2
         p%fine(s)%distribution=dist_normal
         allocate(p%fine(s)%gamma(2,2),p%fine(s)%mu(2),p%fine(s)%sigma(2),p%fine(s)%df(2))
         p%fine(s)%sigma=0.4_dp;p%fine(s)%df=10.0_dp
      end do
      p%fine(1)%gamma=reshape([0.85_dp,0.10_dp,0.15_dp,0.90_dp],[2,2])
      p%fine(2)%gamma=reshape([0.80_dp,0.15_dp,0.20_dp,0.85_dp],[2,2])
      p%fine(1)%mu=[-2.0_dp,0.0_dp];p%fine(2)%mu=[0.0_dp,2.0_dp]
   end subroutine
   subroutine assert_close(x,y,tol,label)
      real(dp),intent(in)::x,y,tol;character(len=*),intent(in)::label
      if(abs(x-y)>tol)then;write(*,*)trim(label),x,y,abs(x-y);error stop 1;end if
   end subroutine
end program test_hhmm
