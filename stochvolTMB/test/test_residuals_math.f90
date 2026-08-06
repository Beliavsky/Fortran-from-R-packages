program test_residuals_math
  use stochvoltmb
  implicit none
  type(sv_rng_state) :: rng
  type(sv_parameters) :: p
  type(sv_simulation) :: sim
  type(sv_fit_result) :: fit
  real(dp), allocatable :: r(:)
  real(dp) :: m,sd,x
  integer :: stat

  call check(abs(normal_quantile(normal_cdf(1.25_dp))-1.25_dp)<1.0e-8_dp,'normal inverse')
  call check(abs(student_t_cdf(0.0_dp,5.0_dp)-0.5_dp)<1.0e-14_dp,'t symmetry')
  call check(abs(skew_normal_cdf(0.0_dp,0.0_dp,1.0_dp,0.0_dp)-0.5_dp)<1.0e-12_dp,'SN normal limit')

  p%sigma_y=0.4_dp; p%sigma_h=0.3_dp; p%phi=0.85_dp
  call rng%seed(564738_8)
  call sim_sv(p,4000,sv_gaussian,rng,sim)
  fit%model=sv_gaussian; fit%params=p
  allocate(fit%h(size(sim%h)),r(size(sim%h)))
  fit%h=sim%h
  call standardized_residuals(sim%y,fit,r,stat)
  m=sum(r)/real(size(r),dp)
  sd=sqrt(sum((r-m)**2)/real(size(r)-1,dp))
  call check(stat==sv_ok .and. abs(m)<0.05_dp,'residual mean')
  call check(abs(sd-1.0_dp)<0.05_dp,'residual standard deviation')
  x=quantile_type7([1.0_dp,2.0_dp,3.0_dp,4.0_dp],0.5_dp)
  call check(abs(x-2.5_dp)<1.0e-14_dp,'R type-7 quantile')
  print '(a)','test_residuals_math: PASS'
contains
  subroutine check(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//message
      error stop 1
    end if
  end subroutine check
end program test_residuals_math
