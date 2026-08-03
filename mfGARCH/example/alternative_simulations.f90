program alternative_simulations
  use, intrinsic :: iso_fortran_env, only : int64
  use mfgarch
  implicit none
  type(mfgarch_model) :: model
  type(mfgarch_simulation) :: rv_simulation, diffusion_simulation
  integer :: status

  model%asymmetric = .false.
  model%gamma = 0.0_dp
  model%k = 3
  model%alpha = 0.06_dp
  model%beta = 0.90_dp
  model%m = -0.2_dp
  model%theta = 0.001_dp
  model%w1 = 1.0_dp
  model%w2 = 3.0_dp

  call simulate_mfgarch_rv_dependent(60, model, 5, 48, .false., rv_simulation, status, &
    seed=1234_int64)
  if (status /= mfgarch_success) error stop 'RV-dependent simulation failed'
  write(*,'(a,es14.6)') 'RV-dependent mean tau: ', &
    sum(rv_simulation%tau)/real(size(rv_simulation%tau),dp)

  call simulate_mfgarch_diffusion(40, model, 0.8_dp, 0.1_dp, 5, 12, &
    diffusion_simulation, status, seed=5678_int64)
  if (status /= mfgarch_success) error stop 'diffusion simulation failed'
  write(*,'(a,es14.6)') 'Diffusion mean realized variance: ', &
    sum(diffusion_simulation%realized_variance)/real(size(diffusion_simulation%returns),dp)
end program alternative_simulations
