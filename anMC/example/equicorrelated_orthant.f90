program equicorrelated_orthant
  use anmc
  implicit none

  integer, parameter :: n = 12
  real(dp) :: mu(n), sigma(n,n), design(n,1), exact_probability
  integer :: i
  type(probability_estimate) :: estimate
  type(simulation_control) :: control

  mu = 0.0_dp
  sigma = 0.5_dp
  do i = 1, n
    sigma(i,i) = 1.0_dp
    design(i,1) = real(i-1,dp) / real(n-1,dp)
  end do

  ! For this equicorrelation-0.5 example at threshold zero,
  ! P(max X > 0) = n/(n+1).
  exact_probability = real(n,dp) / real(n+1,dp)

  control%max_outer = 10000
  control%max_inner = 100
  call seed_fortran_rng(1234)
  estimate = proba_max(0.05_dp, 0.0_dp, mu, sigma, design, q=4, method=0, &
                       algo='ANMC', sim_control=control, &
                       prob_control=genz_bretz(maxpts=100000,abseps=1.0e-5_dp))

  if (.not. estimate%ok) then
    write(*,'(a)') 'anMC failed: '//trim(estimate%message)
    stop 1
  end if

  write(*,'(a,f10.6)') 'estimated P(max X > 0): ', estimate%probability
  write(*,'(a,f10.6)') 'exact probability:       ', exact_probability
  write(*,'(a,i0)')    'active dimensions q:     ', estimate%q
end program equicorrelated_orthant
