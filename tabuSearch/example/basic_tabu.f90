program basic_tabu
   use tabu_search, only : dp, tabu_control, tabu_result, run_tabu_search
   implicit none

   type(tabu_control) :: control
   type(tabu_result) :: result
   integer, allocatable :: best(:)
   integer :: initial(10)

   initial = 0
   control%iters = 30
   control%neigh = 10
   control%list_size = 3
   control%n_restarts = 3
   control%repeat_all = 1
   control%seed = 20260813

   call run_tabu_search(10, objective, result, control, initial)
   best = result%best_configuration()

   write(*,'(a,f10.4)') "best objective: ", result%best_value()
   write(*,'(a,*(i0,1x))') "best configuration: ", best
   write(*,'(a,i0)') "stored iterations: ", result%n_records

contains

   function objective(config) result(value)
      integer, intent(in) :: config(:)
      real(dp) :: value
      real(dp), parameter :: score(10) = [2.0_dp, -1.0_dp, 3.0_dp, 0.5_dp, -2.0_dp, &
         4.0_dp, 1.0_dp, -0.5_dp, 2.5_dp, 1.5_dp]

      value = dot_product(score, real(config, dp))
   end function objective

end program basic_tabu
