program test_reference
   use tabu_search, only : dp, tabu_control, tabu_result, run_tabu_search
   implicit none

   type(tabu_control) :: control
   type(tabu_result) :: result
   integer :: initial(6)
   integer :: expected(8,6)
   real(dp), parameter :: tol = 1.0e-12_dp

   initial = 0
   expected = reshape([ &
      0,0,0,0,0,0, &
      0,0,0,0,0,1, &
      0,0,0,0,1,1, &
      0,0,0,1,1,1, &
      0,0,1,1,1,1, &
      0,1,1,1,1,1, &
      1,1,1,1,1,1, &
      0,1,1,1,1,1  &
      ], [8,6], order=[2,1])

   control%iters = 8
   control%neigh = 6
   control%list_size = 2
   control%n_restarts = 0
   control%repeat_all = 1
   control%seed = 17

   call run_tabu_search(6, weighted_objective, result, control, initial)

   call check(result%n_records == 16, "expected preliminary + diversification histories")
   call check(all(result%config_keep(1:8,:) == expected), "preliminary configuration sequence")
   call check(abs(result%best_value() - 63.0_dp) < tol, "best objective")
   call check(all(result%best_configuration() == 1), "best configuration")

   print '(a)', "test_reference: PASS"

contains

   function weighted_objective(config) result(value)
      integer, intent(in) :: config(:)
      real(dp) :: value
      integer :: i

      value = 0.0_dp
      do i = 1, size(config)
         value = value + real(2**(i-1), dp) * real(config(i), dp)
      end do
   end function weighted_objective

   subroutine check(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         write(*,'(a)') "FAIL: " // trim(message)
         error stop 1
      end if
   end subroutine check

end program test_reference
