program multi_mcga_example
  use mcga, only : dp, multi_mcga_result, multi_mcga_optimize
  implicit none
  type(multi_mcga_result) :: res
  integer :: i

  call multi_mcga_optimize(100, 1, 2, -3.0_dp, 3.0_dp, objectives, res, &
                           crossprob=1.0_dp, mutateprob=0.01_dp, elitism=2, maxiter=100, seed=2026)
  print '(a)', 'top five solutions by upstream multi-MCGA rank score:'
  do i = 1, 5
    print '(f12.6,2x,2es14.5,2x,f8.1)', res%population(i,1), res%costs(i,:), res%ranks(i)
  end do
contains
  subroutine objectives(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    f(1) = (x(1)+1.0_dp)**2
    f(2) = (x(1)-1.0_dp)**2
  end subroutine objectives
end program multi_mcga_example
