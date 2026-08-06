program test_model_fits
  use stochvoltmb
  implicit none
  type(sv_rng_state) :: rng
  type(sv_parameters) :: p
  type(sv_simulation) :: sim
  type(sv_fit_result) :: fit
  type(sv_control) :: ctl
  integer :: model

  p%sigma_y=0.25_dp; p%sigma_h=0.25_dp; p%phi=0.88_dp
  p%df=6.0_dp; p%alpha=-1.0_dp; p%rho=-0.5_dp
  ctl%max_outer_iter=8; ctl%max_inner_iter=45; ctl%compute_covariance=.false.
  do model=sv_student_t,sv_skew_gaussian_leverage
    call rng%seed(30000_8+int(model,8))
    call sim_sv(p,45,model,rng,sim)
    call estimate_parameters(sim%y,model,fit,control=ctl)
    call check(finite_real(fit%objective) .and. fit%objective<1.0e90_dp,'finite fitted objective')
    call check(fit%params%sigma_y>0.0_dp .and. fit%params%sigma_h>0.0_dp,'positive fitted scales')
    call check(abs(fit%params%phi)<1.0_dp,'stationary fitted phi')
    if (model==sv_student_t) call check(fit%params%df>2.0_dp,'valid fitted df')
    if (model==sv_leverage .or. model==sv_skew_gaussian_leverage) then
      call check(abs(fit%params%rho)<1.0_dp,'valid fitted rho')
    end if
  end do
  print '(a)','test_model_fits: PASS'
contains
  subroutine check(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//message
      error stop 1
    end if
  end subroutine check
end program test_model_fits
