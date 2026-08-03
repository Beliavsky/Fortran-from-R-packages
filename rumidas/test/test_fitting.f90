program test_fitting
  use rumidas
  implicit none
  integer,parameter::n=80
  real(dp)::x(n),epsilon(n),true_p(2),mu(n),start(2)
  type(mem_spec)::spec
  type(rumidas_fit_result)::fit
  type(rumidas_fit_control)::control
  real(dp),allocatable::ll0(:)
  integer::i,status

  true_p=[0.20_dp,0.65_dp]
  do i=1,n
    epsilon(i)=0.75_dp+0.5_dp*(0.5_dp+0.5_dp*sin(1.713_dp*real(i,dp)))
  end do
  epsilon=epsilon/(sum(epsilon)/real(n,dp))
  mu(1)=1.0_dp
  x(1)=mu(1)*epsilon(1)
  do i=2,n
    mu(i)=(1.0_dp-sum(true_p))*1.0_dp+true_p(1)*x(i-1)+true_p(2)*mu(i-1)
    x(i)=mu(i)*epsilon(i)
  end do
  spec=mem_spec(RUMIDAS_MEM,0,.false.)
  start=[0.10_dp,0.50_dp]
  call mem_loglik_no_skew(start,x,ll0,status)
  call check(status==0,'initial likelihood')
  control%random_starts=1
  control%max_iterations=250
  control%method='bfgs'
  call umemfit(spec,x,fit,status,start=start,control=control)
  call check(allocated(fit%coefficients),'fit coefficients')
  call check(fit%loglik>=sum(ll0)-1.0e-5_dp,'likelihood improvement')
  call check(fit%coefficients(1)>0.0_dp.and.fit%coefficients(2)>0.0_dp,'positive coefficients')
  call check(sum(fit%coefficients)<1.0_dp,'stationarity')
  call check(allocated(fit%conditional),'fitted prediction')
  call check(all(fit%conditional>0.0_dp),'positive prediction')

  print '(a)', 'test_fitting: PASS'
contains
  subroutine check(condition,message)
    logical,intent(in)::condition
    character(len=*),intent(in)::message
    if(.not.condition) error stop message
  end subroutine check
end program test_fitting
