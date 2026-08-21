program test_tcov_parity
  use flexsurv_kinds, only : dp
  use flexsurv_fit, only : flexsurv_spec, initialize_spec
  use flexsurv_distributions, only : dist_exponential
  use flexsurv_multistate, only : flexsurv_transition, msm_parametric
  use flexsurv_multistate_uncertainty, only : transition_parameters_at_time, &
    transition_spline_model_at_time
  use flexsurv_spline, only : survspline_model,spline_time_identity
  implicit none
  type(flexsurv_spec)::sp
  type(flexsurv_transition)::tr
  real(dp),allocatable::p0(:),p5(:)
  real(dp)::theta(2),want0,want5
  type(survspline_model)::sm
  call initialize_spec(sp,dist_exponential,1,[0.1_dp])
  deallocate(sp%reg(1)%x);allocate(sp%reg(1)%x(1,1));sp%reg(1)%x(1,1)=10.0_dp
  theta=[log(0.1_dp),0.2_dp]
  tr%from=1;tr%to=2;tr%row=1;tr%model_kind=msm_parametric;tr%spec=sp;tr%theta=theta
  tr%tcov_parameter=[1];tr%tcov_column=[1];tr%tcov_rate=[1.0_dp]
  p0=transition_parameters_at_time(tr,0.0_dp);p5=transition_parameters_at_time(tr,5.0_dp)
  want0=exp(theta(1)+10.0_dp*theta(2));want5=exp(theta(1)+15.0_dp*theta(2))
  if(abs(p0(1)-want0)>1.0e-12_dp)error stop 'tcov t0'
  if(abs(p5(1)-want5)>1.0e-12_dp)error stop 'tcov t5'
  if(p5(1)<=p0(1))error stop 'tcov update direction'
  tr%model_kind=2;allocate(tr%spline%model%knots(2),tr%spline%model%gamma(2))
  tr%spline%model%knots=[0.0_dp,1.0_dp];tr%spline%model%gamma=[0.1_dp,1.0_dp]
  tr%spline%model%timescale=spline_time_identity;tr%spline_gamma_tcov_rate=[0.2_dp,-0.1_dp]
  sm=transition_spline_model_at_time(tr,3.0_dp)
  if(maxval(abs(sm%gamma-[0.7_dp,0.7_dp]))>1.0e-12_dp)error stop 'spline tcov gamma'
  print *,'test_tcov_parity: PASS'
end program
