program test_multivariate
  use robstattm, only : dp, covariance_result, projection_result, cov_classic, fast_mve, &
    cov_rob_mm, cov_rob_rocke, init_pp
  implicit none
  integer, parameter :: n = 54, p = 3
  real(dp) :: x(n, p), t
  type(covariance_result) :: classic, mve, mm, rocke
  type(projection_result) :: projection
  integer :: i

  do i = 1, n
    t = real(i, dp)
    x(i, 1) = sin(0.37_dp * t) + 0.15_dp * cos(1.1_dp * t)
    x(i, 2) = 0.65_dp * x(i, 1) + 0.55_dp * cos(0.23_dp * t)
    x(i, 3) = -0.35_dp * x(i, 1) + 0.4_dp * x(i, 2) + 0.3_dp * sin(0.71_dp * t)
  end do
  x(1, :) = [10.0_dp, -8.0_dp, 7.0_dp]
  x(2, :) = [-9.0_dp, 9.0_dp, -6.0_dp]
  x(3, :) = [8.0_dp, 7.0_dp, 9.0_dp]

  call cov_classic(x, classic, correlation=.true.)
  call assert_covariance(classic, p, n, 'classic')
  call fast_mve(x, mve, nsamp=120, seed=4101)
  call assert_covariance(mve, p, n, 'MVE')
  call init_pp(x, projection, minimum_directions=180, seed=4102)
  call assert_true(allocated(projection%covariance), 'projection covariance')
  call cov_rob_mm(x, mm, max_iter=40, seed=4103)
  call assert_covariance(mm, p, n, 'MM-SHR')
  call cov_rob_rocke(x, rocke, max_iter=35, seed=4104)
  call assert_covariance(rocke, p, n, 'Rocke')
  call assert_true(sum(diagonal(mm%covariance)) < sum(diagonal(classic%covariance)), &
    'robust scatter resists contamination')
  print '(a)', 'test_multivariate: PASS'
contains
  subroutine assert_covariance(fit, variables, observations, name)
    type(covariance_result), intent(in) :: fit
    integer, intent(in) :: variables, observations
    character(len=*), intent(in) :: name
    call assert_true(allocated(fit%center), trim(name) // ' center')
    call assert_true(size(fit%center) == variables, trim(name) // ' center size')
    call assert_true(size(fit%covariance, 1) == variables, trim(name) // ' covariance size')
    call assert_true(all(diagonal(fit%covariance) > 0.0_dp), trim(name) // ' positive diagonal')
    if (allocated(fit%distances)) call assert_true(size(fit%distances) == observations, &
      trim(name) // ' distance size')
  end subroutine assert_covariance

  pure function diagonal(a) result(value)
    real(dp), intent(in) :: a(:, :)
    real(dp) :: value(min(size(a, 1), size(a, 2)))
    integer :: j
    do j = 1, size(value)
      value(j) = a(j, j)
    end do
  end function diagonal

  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') 'FAIL: ' // trim(message)
      error stop 1
    end if
  end subroutine assert_true
end program test_multivariate
