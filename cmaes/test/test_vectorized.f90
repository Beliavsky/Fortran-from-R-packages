program test_vectorized
  use cmaes, only : dp, cma_control, cma_result, cma_es
  implicit none
  type(cma_control) :: ctrl
  type(cma_result) :: a, b
  real(dp) :: par(4), lo(4), hi(4)

  par = [3.0_dp, -2.0_dp, 1.0_dp, 4.0_dp]
  lo = -6.0_dp
  hi = 6.0_dp
  ctrl%seed = 314159
  ctrl%maxit = 80
  a = cma_es(par, sphere, lo, hi, ctrl)
  ctrl%vectorized = .true.
  b = cma_es(par, sphere, lo, hi, ctrl, sphere_vec)
  if (maxval(abs(a%par - b%par)) > 1.0e-13_dp) error stop "vectorized trajectory mismatch"
  if (abs(a%value - b%value) > 1.0e-14_dp) error stop "vectorized objective mismatch"
  print '(a,es12.4)', 'vectorized value: ', b%value
contains
  function sphere(x) result(v)
    real(dp), intent(in) :: x(:)
    real(dp) :: v
    v = dot_product(x, x)
  end function sphere

  subroutine sphere_vec(x, v)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(out) :: v(:)
    integer :: j
    do j = 1, size(x, 2)
      v(j) = dot_product(x(:, j), x(:, j))
    end do
  end subroutine sphere_vec
end program test_vectorized
