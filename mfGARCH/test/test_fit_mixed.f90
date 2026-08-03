program test_fit_mixed
  use, intrinsic :: iso_fortran_env, only : int64
  use mfgarch
  implicit none
  type(mfgarch_model) :: truth, start
  type(mfgarch_simulation) :: simulation
  type(mfgarch_fit_control) :: control
  type(mfgarch_fit_result) :: result
  real(dp), allocatable :: covariate_low(:)
  real(dp) :: start_llh
  integer :: status, p, nlow

  truth%asymmetric = .false.
  truth%gamma = 0.0_dp
  truth%k = 3
  truth%alpha = 0.06_dp
  truth%beta = 0.90_dp
  truth%m = 0.0_dp
  truth%theta = 0.25_dp
  truth%w1 = 1.0_dp
  truth%w2 = 3.0_dp
  call simulate_mfgarch(180, truth, 0.75_dp, 0.30_dp, 5, 48, simulation, status, &
    seed=13579_int64)
  call assert_true(status == mfgarch_success, 'mixed simulation status')

  nlow = maxval(simulation%low_frequency_period)
  allocate(covariate_low(nlow))
  do p = 1, nlow
    covariate_low(p) = simulation%covariate((p-1)*5+1)
  end do

  start = truth
  start%mu = 0.10_dp
  start%alpha = 0.03_dp
  start%beta = 0.82_dp
  start%m = -0.2_dp
  start%theta = 0.05_dp
  start%w2 = 2.0_dp
  start_llh = log_likelihood(start, simulation%returns, simulation%low_frequency_period, &
    status, covariate_low)
  control%max_iterations = 220
  control%max_function_evaluations = 20000
  control%multi_start = .false.
  control%compute_inference = .false.
  call fit_mfgarch(simulation%returns, simulation%low_frequency_period, start, result, status, &
    covariate=covariate_low, control=control)
  call assert_true(status == mfgarch_success .or. status == mfgarch_not_converged, 'mixed fit status')
  call assert_true(result%log_likelihood > start_llh, 'mixed likelihood improved')
  call assert_true(result%model%w2 >= 1.0_dp, 'weight shape constraint')
  call assert_true(result%tau_forecast > 0.0_dp, 'positive tau forecast')
  call assert_true(result%variance_ratio >= 0.0_dp, 'nonnegative variance ratio')

  print '(a)', 'test_fit_mixed: PASS'

contains

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//trim(label)
      error stop 1
    end if
  end subroutine assert_true

end program test_fit_mixed
