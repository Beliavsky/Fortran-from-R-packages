program test_negative_neigh
   use tabu_search, only : dp, tabu_control, tabu_result, run_tabu_search
   implicit none

   type(tabu_control) :: control
   type(tabu_result) :: result
   integer :: initial(10), i
   real(dp) :: expected

   initial = 0
   control%iters = 40
   control%neigh = 3
   control%list_size = 2
   control%n_restarts = 2
   control%repeat_all = 1
   control%seed = 987654321

   call run_tabu_search(10, negative_objective, result, control, initial)

   do i = 1, result%n_records
      expected = negative_objective(result%config_keep(i,:))
      call check(abs(result%utility_keep(i) - expected) < 1.0e-12_dp, &
         "stored utility must correspond to the stored configuration")
   end do
   call check(all(result%config_keep == 0 .or. result%config_keep == 1), "binary histories")
   call check(result%best_value() <= 0.0_dp, "negative objective remains negative or zero")

   print '(a)', "test_negative_neigh: PASS"

contains

   function negative_objective(config) result(value)
      integer, intent(in) :: config(:)
      real(dp) :: value
      integer, parameter :: target(10) = [1,0,1,1,0,1,0,1,0,1]

      value = -real(count(config /= target), dp)
   end function negative_objective

   subroutine check(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         write(*,'(a)') "FAIL: " // trim(message)
         error stop 1
      end if
   end subroutine check

end program test_negative_neigh
