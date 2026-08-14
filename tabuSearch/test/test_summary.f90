program test_summary
   use tabu_search, only : dp, tabu_control, tabu_result, tabu_summary_result, &
      run_tabu_search, summarize_tabu
   implicit none

   type(tabu_control) :: control
   type(tabu_result) :: result
   type(tabu_summary_result) :: summary
   integer :: initial(7), expected_moves(7), i

   initial = 0
   control%iters = 12
   control%neigh = 7
   control%list_size = 2
   control%n_restarts = 1
   control%repeat_all = 1
   control%seed = 41

   call run_tabu_search(7, objective, result, control, initial)
   call summarize_tabu(result, summary)

   call check(summary%total_iterations == result%n_records, "summary iteration count")
   call check(summary%best_count >= 1, "at least one optimum occurrence")
   call check(summary%unique_configurations >= 1, "unique configuration count")
   call check(size(summary%selected_count) == 7, "selected-count size")
   call check(size(summary%move_frequency) == 7, "move-frequency size")
   expected_moves = 0
   do i = 2, result%n_records
      where (result%config_keep(i,:) /= result%config_keep(i-1,:))
         expected_moves = expected_moves + 1
      end where
   end do
   call check(all(summary%move_frequency == expected_moves), "move-frequency values")
   call check(abs(summary%best_value - result%best_value()) < 1.0e-12_dp, "summary best value")

   print '(a)', "test_summary: PASS"

contains

   function objective(config) result(value)
      integer, intent(in) :: config(:)
      real(dp) :: value

      value = real(sum(config), dp)
   end function objective

   subroutine check(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         write(*,'(a)') "FAIL: " // trim(message)
         error stop 1
      end if
   end subroutine check

end program test_summary
