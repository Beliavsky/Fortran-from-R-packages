module test_optimization_functions
  use multiatsm_kinds, only : dp
  implicit none
contains
  function quadratic(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    value = 0.5_dp * (x(1) - 1.5_dp)**2 + 2.0_dp * (x(2) + 0.75_dp)**2
  end function quadratic

  subroutine vector_test(x, value)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: value(:)
    allocate(value(2))
    value = [x(1)**2 + x(2), x(1) * x(2)]
  end subroutine vector_test
end module test_optimization_functions

program test_optimization
  use multiatsm_kinds, only : dp
  use multiatsm_linalg, only : spectral_radius
  use multiatsm_optimization, only : optimization_result, bfgs_minimize, nelder_mead_minimize, &
    numerical_gradient, numerical_jacobian, stabilize_transition, lower_factor_to_psd, &
    psd_to_lower_parameters, block_diagonal_psd
  use test_optimization_functions, only : quadratic, vector_test
  implicit none
  real(dp) :: x0(2), matrix(2, 2), radius
  real(dp), allocatable :: gradient(:), jacobian(:, :), stable(:, :), covariance(:, :), lower(:, :)
  real(dp), allocatable :: parameters(:), block_covariance(:, :)
  type(optimization_result) :: bfgs, nm
  integer :: info

  x0 = [4.0_dp, -3.0_dp]
  call bfgs_minimize(quadratic, x0, bfgs, tolerance=1.0e-10_dp, max_iterations=300)
  call check(bfgs%converged, 'BFGS convergence')
  call check(maxval(abs(bfgs%x - [1.5_dp, -0.75_dp])) < 1.0e-5_dp, 'BFGS optimum')
  call nelder_mead_minimize(quadratic, x0, nm, tolerance=1.0e-9_dp, max_iterations=1000)
  call check(nm%converged, 'Nelder-Mead convergence')
  call check(maxval(abs(nm%x - [1.5_dp, -0.75_dp])) < 1.0e-4_dp, 'Nelder-Mead optimum')

  call numerical_gradient(quadratic, [2.0_dp, 1.0_dp], gradient)
  call check(maxval(abs(gradient - [0.5_dp, 7.0_dp])) < 1.0e-5_dp, 'numerical gradient')
  call numerical_jacobian(vector_test, [2.0_dp, 3.0_dp], jacobian, info)
  call check(info == 0, 'Jacobian status')
  call check(maxval(abs(jacobian - reshape([4.0_dp, 3.0_dp, 1.0_dp, 2.0_dp], [2, 2]))) < 1.0e-5_dp, &
    'numerical Jacobian')

  matrix = reshape([1.2_dp, 0.0_dp, 0.1_dp, 0.8_dp], [2, 2])
  call stabilize_transition(matrix, stable, info, 0.95_dp)
  call check(info == 0, 'stabilization status')
  radius = spectral_radius(stable, info)
  call check(info == 0 .and. radius <= 0.9500001_dp, 'stationarity restriction')

  parameters = [1.0_dp, 0.2_dp, 0.7_dp]
  call lower_factor_to_psd(parameters, 2, covariance, lower, info)
  call check(info == 0, 'PSD transform status')
  call check(maxval(abs(covariance - matmul(lower, transpose(lower)))) < 1.0e-12_dp, 'PSD factorization')
  call psd_to_lower_parameters(covariance, parameters, info)
  call check(info == 0 .and. size(parameters) == 3, 'PSD auxiliary parameters')
  call lower_factor_to_psd(parameters, 2, block_covariance, lower, info)
  call check(info == 0 .and. maxval(abs(block_covariance - covariance)) < 1.0e-9_dp, &
    'PSD parameter round trip')
  call block_diagonal_psd([1.0_dp, 0.2_dp, 0.7_dp, 0.5_dp], [2, 1], block_covariance, info)
  call check(info == 0, 'block PSD status')
  call check(maxval(abs(block_covariance(1:2, 3))) < 1.0e-12_dp, 'block covariance zeros')
  print '(a)', 'test_optimization: PASS'
contains
  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') 'FAIL: ' // message
      error stop 1
    end if
  end subroutine check
end program test_optimization
