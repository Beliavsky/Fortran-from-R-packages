program fit_simple_garch
  use, intrinsic :: iso_fortran_env, only : int64
  use mfgarch
  implicit none
  type(mfgarch_model) :: truth, start
  type(mfgarch_simulation) :: simulation
  type(mfgarch_fit_control) :: control
  type(mfgarch_fit_result) :: fit
  integer, allocatable :: period(:)
  integer :: status, i

  truth%asymmetric = .false.
  truth%gamma = 0.0_dp
  truth%k = 0
  truth%alpha = 0.07_dp
  truth%beta = 0.90_dp
  truth%m = log(1.2_dp)
  call simulate_mfgarch(300, truth, 0.0_dp, 0.0_dp, 1, 48, simulation, status, &
    seed=112233_int64)
  allocate(period(300))
  period = [(i,i=1,300)]

  start = truth
  start%alpha = 0.03_dp
  start%beta = 0.80_dp
  start%m = 0.0_dp
  control%multi_start = .false.
  call fit_mfgarch(simulation%returns, period, start, fit, status, control=control)
  call print_fit_summary(fit)
end program fit_simple_garch
