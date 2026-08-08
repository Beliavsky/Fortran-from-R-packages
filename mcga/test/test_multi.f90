program test_multi
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mcga, only : dp, multi_mcga_result, multi_mcga_optimize
  implicit none
  type(multi_mcga_result) :: res
  integer :: i

  call multi_mcga_optimize(80, 1, 2, -4.0_dp, 4.0_dp, objectives, res, &
                           crossprob=1.0_dp, mutateprob=0.01_dp, elitism=2, maxiter=80, seed=1357)
  if (size(res%population,1) /= 80 .or. size(res%costs,2) /= 2) error stop "multi result shape failed"
  do i = 2, size(res%ranks)
    if (res%ranks(i) > res%ranks(i-1)) error stop "multi output not rank-sorted"
  end do
  if (.not. all(ieee_is_finite(res%costs))) error stop "multi costs are not finite"
  print *, "test_multi: PASS", res%population(1,1), res%ranks(1)
contains
  subroutine objectives(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    f(1) = (x(1) + 1.0_dp)**2
    f(2) = (x(1) - 1.0_dp)**2
  end subroutine objectives
end program test_multi
