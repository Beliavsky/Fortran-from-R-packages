program demo_mfgarch
  use, intrinsic :: iso_fortran_env, only : int64
  use mfgarch
  implicit none
  type(mfgarch_model) :: data_model, start_model
  type(mfgarch_simulation) :: simulation
  type(mfgarch_fit_control) :: control
  type(mfgarch_fit_result) :: fit
  real(dp), allocatable :: covariate_low(:), forecasts(:)
  integer :: status, p, nlow

  data_model%asymmetric = .true.
  data_model%k = 6
  data_model%alpha = 0.05_dp
  data_model%beta = 0.88_dp
  data_model%gamma = 0.08_dp
  data_model%m = 0.0_dp
  data_model%theta = 0.30_dp
  data_model%w1 = 1.0_dp
  data_model%w2 = 4.0_dp

  call simulate_mfgarch(300, data_model, 0.85_dp, 0.25_dp, 5, 48, simulation, status, &
    seed=20260801_int64)
  if (status /= mfgarch_success) error stop 'simulation failed'

  nlow = maxval(simulation%low_frequency_period)
  allocate(covariate_low(nlow))
  do p = 1, nlow
    covariate_low(p) = simulation%covariate((p-1)*5+1)
  end do

  start_model = data_model
  start_model%mu = 0.0_dp
  start_model%alpha = 0.03_dp
  start_model%beta = 0.85_dp
  start_model%gamma = 0.04_dp
  start_model%theta = 0.10_dp
  start_model%w2 = 3.0_dp
  control%max_iterations = 300
  control%multi_start = .false.
  control%compute_inference = .true.

  call fit_mfgarch(simulation%returns, simulation%low_frequency_period, start_model, fit, &
    status, covariate=covariate_low, control=control)
  call print_fit_summary(fit)

  call predict_variance(fit%model, [1,2,5,10], fit%tau_forecast, &
    simulation%returns(size(simulation%returns)), fit%g(size(fit%g)), &
    fit%tau(size(fit%tau)), forecasts, status)
  if (status == mfgarch_success) then
    write(*,'(/,a)') 'Variance forecasts:'
    write(*,'(*(es14.6,1x))') forecasts
  end if
end program demo_mfgarch
