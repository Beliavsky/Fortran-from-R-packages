program mcga_example
  use mcga, only : dp, mcga_result, mcga_optimize
  implicit none
  type(mcga_result) :: res

  call mcga_optimize(popsize=150, chsize=3, minval=-100.0_dp, maxval=100.0_dp, &
                     objective=objective, result=res, crossprob=1.0_dp, mutateprob=0.01_dp, &
                     elitism=2, maxiter=250, seed=2026)
  print '(a,3f14.6)', 'best parameters: ', res%population(1,:)
  print '(a,es14.6)', 'best cost:       ', res%costs(1)
contains
  function objective(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    f = (x(1)-7.0_dp)**2 + (x(2)-17.0_dp)**2 + (x(3)-27.0_dp)**2
  end function objective
end program mcga_example
