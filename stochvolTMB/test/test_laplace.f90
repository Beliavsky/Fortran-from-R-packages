program test_laplace
  use stochvoltmb
  implicit none
  type(sv_rng_state) :: rng
  type(sv_parameters) :: p
  type(sv_simulation) :: sim
  type(sv_control) :: ctl
  real(dp), allocatable :: theta(:), h(:), d(:), o(:)
  real(dp) :: value, corr
  integer :: model, stat

  p%sigma_y=0.25_dp; p%sigma_h=0.3_dp; p%phi=0.92_dp
  p%df=7.0_dp; p%alpha=-1.5_dp; p%rho=-0.55_dp
  ctl%max_inner_iter=80; ctl%compute_covariance=.false.
  do model=sv_gaussian,sv_skew_gaussian_leverage
    call rng%seed(20000_8+int(model,8))
    call sim_sv(p,100,model,rng,sim)
    allocate(theta(parameter_count(model)),h(100),d(100),o(99))
    call parameters_to_theta(p,model,theta,stat)
    call check(stat==sv_ok,'parameter transform')
    value=laplace_nll(sim%y,theta,model,ctl,h,d,o,stat)
    call check(finite_real(value) .and. value<1.0e90_dp,'finite Laplace objective')
    call check(all(d>0.0_dp),'positive latent Hessian diagonal')
    if (model==sv_gaussian) then
      corr=correlation(h,sim%h)
      call check(corr>0.25_dp,'latent mode tracks simulated volatility')
    end if
    deallocate(theta,h,d,o)
  end do
  print '(a)','test_laplace: PASS'
contains
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
end program test_laplace
