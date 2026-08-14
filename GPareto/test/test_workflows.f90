program test_workflows
  use gpareto, only : dp, gp_model_set, fit_gp_model, trend_const
  use gpareto, only : gpareto_result, gparetoptim, optim_control
  use gpareto, only : pareto_density_result, pareto_set_density, integration_design_optim
  implicit none
  type(gp_model_set)::ms
  type(gpareto_result)::res
  type(pareto_density_result)::dens
  type(optim_control)::ctl
  real(dp)::x(7,1),y1(7),y2(7),ref(2)
  real(dp),allocatable::ip(:,:),iw(:)
  integer::i
  do i=1,7
    x(i,1)=real(i-1,dp)/6.0_dp
    y1(i)=x(i,1)**2
    y2(i)=(1.0_dp-x(i,1))**2
  end do
  allocate(ms%model(2))
  call fit_gp_model(ms%model(1),x,y1,covtype='gauss',trend_kind=trend_const,nugget=1.0e-8_dp)
  call fit_gp_model(ms%model(2),x,y2,covtype='gauss',trend_kind=trend_const,nugget=1.0e-8_dp)
  ctl%population=10
  ctl%generations=5
  ctl%seed=19
  ref=[1.2_dp,1.2_dp]
  call gparetoptim(ms,obj,1,[0.0_dp],[1.0_dp],'SMS',res,ref=ref, &
    cov_reestimate=.false.,control=ctl)
  if(size(res%x_new,1)/=1 .or. ms%model(1)%km%n/=8) error stop 'gparetoptim failed'
  if(any(res%x_new<0.0_dp).or.any(res%x_new>1.0_dp)) error stop 'gparetoptim bounds failed'
  call integration_design_optim([0.0_dp],[1.0_dp],ip,iw,npoints=12,distribution='halton')
  if(size(ip,1)/=12.or.size(iw)/=0) error stop 'integration design failed'
  call integration_design_optim([0.0_dp],[1.0_dp],ip,iw,npoints=8,distribution='SUR', &
    models=ms,front=res%front,seed=31_8)
  if(size(ip,1)/=8.or.size(iw)/=8.or.any(iw<=0.0_dp)) error stop 'SUR integration design failed'
  call pareto_set_density(ms,[0.0_dp],[1.0_dp],dens,nsim=4,npoints=20,seed=33_8)
  if(size(dens%cps,1)==0 .or. any(dens%density<0.0_dp)) error stop 'ParetoSetDensity failed'
  print *, 'test_workflows PASS'
contains
  subroutine obj(xx,yy)
    real(dp),intent(in)::xx(:)
    real(dp),intent(out)::yy(:)
    yy(1)=xx(1)**2
    yy(2)=(1.0_dp-xx(1))**2
  end subroutine obj
end program test_workflows
