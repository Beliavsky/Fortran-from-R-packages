program test_garch_fitting
  use rumidas
  implicit none
  integer,parameter::n=100,k=2
  real(dp)::ret(n),z(n),g(n),tau,mv(k+1,n),ptrue(5),start(5)
  type(garch_midas_spec)::spec
  type(rumidas_fit_result)::fit
  type(rumidas_fit_control)::control
  real(dp),allocatable::ll0(:)
  integer::i,status

  ptrue=[0.08_dp,0.84_dp,log(1.0e-4_dp),0.0_dp,2.0_dp]
  tau=exp(ptrue(3));mv=0.2_dp
  do i=1,n
    z(i)=sin(1.37_dp*real(i,dp))+0.4_dp*cos(0.53_dp*real(i,dp))
  end do
  z=z/sqrt(sum(z*z)/real(n,dp))
  g(1)=1.0_dp;ret(1)=sqrt(tau*g(1))*z(1)
  do i=2,n
    g(i)=1.0_dp-ptrue(1)-ptrue(2)+ptrue(1)*ret(i-1)**2/tau+ptrue(2)*g(i-1)
    ret(i)=sqrt(tau*g(i))*z(i)
  end do
  spec=garch_midas_spec(RUMIDAS_GM,RUMIDAS_NORMAL,RUMIDAS_BETA_LAG,k,0,.false.)
  start=[0.05_dp,0.75_dp,log(sum(ret*ret)/real(n,dp)),0.0_dp,2.0_dp]
  call gm_loglik_no_skew(start,ret,mv,k,ll0,status)
  call check(status==0,'initial evaluation')
  control%random_starts=1;control%max_iterations=300;control%method='bfgs'
  call ugmfit(spec,ret,mv,fit,status,start=start,control=control)
  call check(allocated(fit%coefficients),'coefficients')
  call check(fit%loglik>=sum(ll0)-1.0e-4_dp,'likelihood improvement')
  call check(fit%coefficients(1)>0.0_dp.and.fit%coefficients(2)>0.0_dp,'positive')
  call check(fit%coefficients(1)+fit%coefficients(2)<1.0_dp,'stationary')
  call check(all(fit%conditional>0.0_dp),'positive variance')
  print '(a)', 'test_garch_fitting: PASS'
contains
  subroutine check(condition,message)
    logical,intent(in)::condition
    character(len=*),intent(in)::message
    if(.not.condition) error stop message
  end subroutine check
end program test_garch_fitting
