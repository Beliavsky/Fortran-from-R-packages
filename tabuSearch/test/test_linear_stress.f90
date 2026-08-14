program test_linear_stress
   use tabu_search, only : dp, i8, tabu_control, tabu_result, run_tabu_search
   use tabu_search_rng, only : tabu_rng
   implicit none

   integer, parameter :: n = 20, ncases = 100
   type(tabu_control) :: control
   type(tabu_result) :: result
   type(tabu_rng) :: data_rng
   real(dp) :: weights(n), expected, u
   integer :: initial(n), optimum(n), case_index, i

   call data_rng%seed(314159_i8)
   control%iters = n + 3
   control%neigh = n
   control%list_size = 4
   control%n_restarts = 0
   control%repeat_all = 1

   do case_index = 1, ncases
      do i = 1, n
         u = data_rng%uniform()
         weights(i) = 4.0_dp * u - 2.0_dp
         if (abs(weights(i)) < 0.05_dp) weights(i) = sign(0.05_dp, weights(i) + 1.0e-12_dp)
         initial(i) = data_rng%randint(0, 1)
      end do
      optimum = 0
      where (weights > 0.0_dp) optimum = 1
      expected = dot_product(weights, real(optimum, dp))
      control%seed = int(1000 + case_index, i8)

      call run_tabu_search(n, objective, result, control, initial)
      call check(abs(result%best_value() - expected) < 1.0e-10_dp, "linear optimum value")
      call check(all(result%best_configuration() == optimum), "linear optimum configuration")
   end do

   print '(a,i0,a)', "test_linear_stress: PASS (", ncases, " cases)"

contains

   function objective(config) result(value)
      integer, intent(in) :: config(:)
      real(dp) :: value

      value = dot_product(weights, real(config, dp))
   end function objective

   subroutine check(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         write(*,'(a)') "FAIL: " // trim(message)
         error stop 1
      end if
   end subroutine check

end program test_linear_stress
