program test_fit
  use stochvoltmb
  implicit none
  type(sv_rng_state) :: rng
  type(sv_parameters) :: p,start
  type(sv_simulation) :: sim
  type(sv_fit_result) :: fit
  type(sv_control) :: ctl
  real(dp), allocatable :: theta0(:)
  real(dp) :: initial_objective
  integer :: stat

  p%sigma_y=0.3_dp; p%sigma_h=0.28_dp; p%phi=0.93_dp
  call rng%seed(777777_8)
  call sim_sv(p,160,sv_gaussian,rng,sim)
  start%sigma_y=sample_sd(sim%y); start%sigma_h=0.25_dp; start%phi=0.95_dp
  allocate(theta0(3))
  call parameters_to_theta(start,sv_gaussian,theta0,stat)
  ctl%max_outer_iter=140; ctl%max_inner_iter=70
  ctl%compute_covariance=.false.; ctl%outer_tolerance=2.0e-5_dp
  initial_objective=get_nll(sim%y,theta0,sv_gaussian,ctl)
  call estimate_parameters(sim%y,sv_gaussian,fit,start=start,control=ctl)
  call check(fit%status==sv_ok .or. fit%status==sv_no_convergence,'fit status')
  call check(fit%objective<=initial_objective+1.0e-5_dp,'optimizer improves objective')
  call check(fit%params%sigma_y>0.0_dp .and. fit%params%sigma_h>0.0_dp,'positive scales')
  call check(abs(fit%params%phi)<1.0_dp,'stationary persistence')
  call check(size(fit%h)==160 .and. all(fit%h_se>=0.0_dp),'latent output')
  call check(correlation(fit%h,sim%h)>0.25_dp,'fitted latent volatility')
  print '(a)','test_fit: PASS'
contains
  real(dp) function sample_sd(x) result(s)
    real(dp), intent(in) :: x(:)
    real(dp) :: m
    m=sum(x)/real(size(x),dp)
    s=sqrt(sum((x-m)**2)/real(size(x)-1,dp))
  end function sample_sd
  real(dp) function correlation(x,y) result(r)
    real(dp), intent(in) :: x(:),y(:)
    real(dp) :: xm,ym
    xm=sum(x)/real(size(x),dp); ym=sum(y)/real(size(y),dp)
    r=sum((x-xm)*(y-ym))/sqrt(sum((x-xm)**2)*sum((y-ym)**2))
  end function correlation
  subroutine check(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//message
      error stop 1
    end if
  end subroutine check
end program test_fit
