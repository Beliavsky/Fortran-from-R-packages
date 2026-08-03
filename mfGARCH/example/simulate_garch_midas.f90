program simulate_garch_midas
  use, intrinsic :: iso_fortran_env, only : int64
  use mfgarch
  implicit none
  type(mfgarch_model) :: model
  type(mfgarch_simulation) :: simulation
  integer :: status

  model%asymmetric = .true.
  model%k = 12
  model%alpha = 0.06_dp
  model%beta = 0.88_dp
  model%gamma = 0.06_dp
  model%m = 0.0_dp
  model%theta = 0.20_dp
  model%w1 = 1.0_dp
  model%w2 = 3.0_dp

  call simulate_mfgarch(250, model, 0.95_dp, 0.10_dp, 5, 48, simulation, status, &
    seed=24680_int64, student_t_df=8.0_dp)
  if (status /= mfgarch_success) error stop 'simulation failed'
  call write_simulation_csv('mfgarch_simulation.csv', simulation, status)
  write(*,'(a,f12.6)') 'Mean daily return: ', sum(simulation%returns)/real(size(simulation%returns),dp)
  write(*,'(a,f12.6)') 'Mean conditional variance: ', &
    sum(simulation%tau*simulation%g)/real(size(simulation%returns),dp)
  write(*,'(a)') 'Wrote mfgarch_simulation.csv'
end program simulate_garch_midas
