program vectorized_example
  use cmaes, only : dp, cma_control, cma_result, cma_es
  implicit none
  type(cma_control) :: control
  type(cma_result) :: result
  real(dp) :: par(3), lower(3), upper(3)

  par = [2.0_dp, -3.0_dp, 4.0_dp]
  lower = -5.0_dp
  upper = 5.0_dp
  control%seed = 8128
  control%maxit = 300
  control%vectorized = .true.
  control%stopfitness = 1.0e-10_dp

  result = cma_es(par, scalar_sphere, lower, upper, control, vector_sphere)
  print '(a,es16.8)', 'value = ', result%value
  print '(a,*(f12.6,1x))', 'par   = ', result%par
contains
  function scalar_sphere(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    value = dot_product(x, x)
  end function scalar_sphere

  subroutine vector_sphere(x, value)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(out) :: value(:)
    integer :: j
    do j = 1, size(x, 2)
      value(j) = dot_product(x(:, j), x(:, j))
    end do
  end subroutine vector_sphere
end program vectorized_example
