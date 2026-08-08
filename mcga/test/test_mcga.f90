program test_mcga
  use mcga, only : dp, mcga_result, mcga_optimize
  implicit none
  type(mcga_result) :: res

  call mcga_optimize(120, 2, -10.0_dp, 10.0_dp, objective, res, &
                     crossprob=1.0_dp, mutateprob=0.01_dp, elitism=2, maxiter=180, seed=2468)
  if (res%costs(1) > 1.0e-2_dp) then
    print *, "best = ", res%population(1,:), " cost = ", res%costs(1)
    error stop "mcga quadratic optimization failed"
  end if
  if (any(res%costs(2:) < res%costs(:size(res%costs)-1))) error stop "mcga output not sorted"
  print *, "test_mcga: PASS", res%population(1,:), res%costs(1)
contains
  function objective(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    f = (x(1) - 2.5_dp)**2 + (x(2) + 1.25_dp)**2
  end function objective
end program test_mcga
