program test_tail_ratios
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use rpeif, only : dp, rpeif_options, nuisance_parameters, nuisance_parameters_fn, &
    influence_from_data, influence_from_nuisance, rpeif_success, rpeif_numerical_failure
  implicit none
  real(dp) :: returns(12), x_eval(7)
  real(dp), allocatable :: values(:), source_values(:), corrected_values(:)
  type(rpeif_options) :: opts
  type(nuisance_parameters) :: pars
  character(len=16), parameter :: estimators(6) = [character(len=16) :: &
    'var', 'es', 'esratio', 'varratio', 'rachevratio', 'robmean']
  integer :: i, status

  returns = [-0.12_dp, -0.08_dp, -0.04_dp, -0.02_dp, -0.005_dp, 0.01_dp, &
    0.02_dp, 0.035_dp, 0.05_dp, 0.08_dp, 0.11_dp, 0.16_dp]
  x_eval = [-0.15_dp, -0.08_dp, -0.02_dp, 0.0_dp, 0.04_dp, 0.10_dp, 0.18_dp]
  opts%alpha = 0.2_dp
  opts%beta = 0.2_dp

  do i = 1, size(estimators)
    call influence_from_data(trim(estimators(i)), x_eval, returns, values, opts, status)
    call assert_true(status == rpeif_success .or. status == rpeif_numerical_failure, trim(estimators(i))//' status')
    call assert_true(all(ieee_is_finite(values)), trim(estimators(i))//' finite')
  end do

  opts%source_compatibility = .true.
  call influence_from_data('varratio', x_eval, returns, source_values, opts, status)
  opts%source_compatibility = .false.
  call influence_from_data('varratio', x_eval, returns, corrected_values, opts, status)
  call assert_true(maxval(abs(source_values - corrected_values)) > 1.0e-6_dp, 'VaR-ratio compatibility switch')

  call nuisance_parameters_fn(pars, mu=0.01_dp, sd=0.05_dp, alpha=0.2_dp, beta=0.2_dp, status=status)
  call influence_from_nuisance('rachevratio', x_eval, pars, values, opts, status)
  call assert_true(status == rpeif_success, 'nuisance Rachev status')
  call assert_true(all(ieee_is_finite(values)), 'nuisance Rachev finite')

  print '(a)', 'test_tail_ratios: PASS'
contains
  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      print '(a)', trim(label)//' failed'
      error stop 1
    end if
  end subroutine assert_true
end program test_tail_ratios
