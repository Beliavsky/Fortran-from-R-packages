program test_diagnostics
  use cmaes, only : dp, cma_control, cma_result, cma_es, extract_population
  implicit none
  type(cma_control) :: ctrl
  type(cma_result) :: res
  real(dp) :: par(3), lo(3), hi(3)
  real(dp), allocatable :: pop(:, :), values(:)

  par = [2.0_dp, -1.0_dp, 3.0_dp]
  lo = -5.0_dp
  hi = 5.0_dp
  ctrl%seed = 911
  ctrl%maxit = 30
  ctrl%diag = .true.
  res = cma_es(par, sphere, lo, hi, ctrl)
  if (.not. allocated(res%sigma_history)) error stop "missing sigma history"
  if (.not. allocated(res%eigen_history)) error stop "missing eigen history"
  if (.not. allocated(res%value_history)) error stop "missing value history"
  if (.not. allocated(res%population_history)) error stop "missing population history"
  if (size(res%sigma_history) /= res%iterations) error stop "sigma history size"
  call extract_population(res, min(2, res%iterations), pop, values)
  if (size(pop, 1) /= 3) error stop "population dimension"
  if (size(values) /= size(pop, 2)) error stop "population values dimension"
  print '(a,i0)', 'diagnostic iterations: ', res%iterations
contains
  function sphere(x) result(v)
    real(dp), intent(in) :: x(:)
    real(dp) :: v
    v = dot_product(x, x)
  end function sphere
end program test_diagnostics
