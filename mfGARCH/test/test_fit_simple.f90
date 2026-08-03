program test_fit_simple
  use, intrinsic :: iso_fortran_env, only : int64
  use mfgarch
  implicit none
  type(mfgarch_model) :: truth, start
  type(mfgarch_simulation) :: simulation
  type(mfgarch_fit_control) :: control
  type(mfgarch_fit_result) :: result
  integer, allocatable :: period(:)
  real(dp) :: start_llh
  integer :: status, i

  truth%asymmetric = .false.
  truth%gamma = 0.0_dp
  truth%k = 0
  truth%alpha = 0.08_dp
  truth%beta = 0.88_dp
  truth%m = log(1.3_dp)
  call simulate_mfgarch(240, truth, 0.0_dp, 0.0_dp, 1, 48, simulation, status, &
    seed=20260801_int64)
  call assert_true(status == mfgarch_success, 'simulation status')

  allocate(period(240))
  period = [(i,i=1,240)]
  start = truth
  start%mu = 0.15_dp
  start%alpha = 0.03_dp
  start%beta = 0.75_dp
  start%m = -0.4_dp
  start_llh = log_likelihood(start, simulation%returns, period, status)
  control%max_iterations = 250
  control%max_function_evaluations = 15000
  control%multi_start = .false.
  control%compute_inference = .true.
  call fit_mfgarch(simulation%returns, period, start, result, status, control=control)
  call assert_true(status == mfgarch_success .or. status == mfgarch_not_converged, 'fit status')
  call assert_true(result%log_likelihood > start_llh, 'likelihood improved')
  call assert_true(result%model%alpha >= 0.0_dp .and. result%model%beta >= 0.0_dp, &
    'nonnegative GARCH parameters')
  call assert_true(result%model%alpha + result%model%beta < 1.0_dp, 'stationarity')
  call assert_true(size(result%standard_error) == 4, 'standard errors available')
  call assert_true(all(result%standard_error >= 0.0_dp), 'standard errors nonnegative')
  call assert_true(all(result%tau > 0.0_dp) .and. all(result%g > 0.0_dp), 'fitted components positive')

  print '(a)', 'test_fit_simple: PASS'

contains

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//trim(label)
      error stop 1
    end if
  end subroutine assert_true

end program test_fit_simple
