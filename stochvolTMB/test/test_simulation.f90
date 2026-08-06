program test_simulation
  use stochvoltmb
  implicit none
  type(sv_rng_state) :: rng
  type(sv_parameters) :: p
  type(sv_simulation) :: sim
  integer :: model
  real(dp) :: m, v

  p%sigma_y=0.3_dp; p%sigma_h=0.25_dp; p%phi=0.9_dp
  p%df=6.0_dp; p%alpha=-2.0_dp; p%rho=-0.6_dp
  do model=sv_gaussian,sv_skew_gaussian_leverage
    call rng%seed(10000_8+int(model,8))
    call sim_sv(p,600,model,rng,sim)
    call check(sim%status==sv_ok,'simulation status')
    call check(size(sim%y)==600 .and. size(sim%h)==600,'simulation dimensions')
    call check(all(finite_real(sim%y)) .and. all(finite_real(sim%h)),'finite simulation')
  end do

  call rng%seed(12345_8)
  p%sigma_h=0.0_dp; p%phi=0.0_dp; p%alpha=4.0_dp
  call sim_sv(p,20000,sv_skew_gaussian,rng,sim)
  m=sum(sim%y)/real(size(sim%y),dp)
  v=sum((sim%y-m)**2)/real(size(sim%y)-1,dp)
  call check(abs(m)<0.015_dp,'standardized skew-normal mean')
  call check(abs(v-p%sigma_y**2)<0.01_dp,'standardized skew-normal variance')

  p%phi=1.1_dp
  call sim_sv(p,100,sv_gaussian,rng,sim)
  call check(sim%status==sv_invalid_argument,'invalid phi')
  print '(a)','test_simulation: PASS'
contains
  subroutine check(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//message
      error stop 1
    end if
  end subroutine check
end program test_simulation
