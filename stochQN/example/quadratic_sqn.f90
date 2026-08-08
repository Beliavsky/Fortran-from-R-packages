program quadratic_sqn
   use stochqn_kinds, only : dp
   use stochqn_guided, only : optimize_sqn, stochqn_run_result_t
   implicit none

   type :: problem_t
      real(dp) :: diagonal(3)
      real(dp) :: rhs(3)
   end type problem_t

   type(problem_t) :: problem
   type(stochqn_run_result_t) :: result
   real(dp) :: x(3)

   problem%diagonal = [1.0_dp, 2.0_dp, 4.0_dp]
   problem%rhs = [1.0_dp, -2.0_dp, 0.5_dp]
   x = [4.0_dp, -3.0_dp, 2.0_dp]

   call optimize_sqn(x, 400, 0.12_dp, gradient, result, hessian_vector, problem, &
                     bfgs_update_frequency=5)

   print '(a,3f14.8)', 'solution: ', x
   print '(a,i0)', 'iterations: ', result%iterations
   print '(a,i0)', 'gradient evaluations: ', result%gradient_evaluations
   print '(a,i0)', 'Hessian-vector evaluations: ', result%hessian_vector_evaluations

contains

   subroutine gradient(x, task, iteration, g, user_data)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: task, iteration
      real(dp), intent(out) :: g(:)
      class(*), intent(inout), optional :: user_data

      select type (data => user_data)
      type is (problem_t)
         g = data%diagonal * x - data%rhs
      class default
         g = 0.0_dp
      end select
      if (task + iteration < -huge(1)) g = 0.0_dp
   end subroutine gradient

   subroutine hessian_vector(x, vector, iteration, product, user_data)
      real(dp), intent(in) :: x(:), vector(:)
      integer, intent(in) :: iteration
      real(dp), intent(out) :: product(:)
      class(*), intent(inout), optional :: user_data

      select type (data => user_data)
      type is (problem_t)
         product = data%diagonal * vector
      class default
         product = 0.0_dp
      end select
      if (sum(x) + real(iteration, dp) < -huge(1.0_dp)) product = 0.0_dp
   end subroutine hessian_vector

end program quadratic_sqn
