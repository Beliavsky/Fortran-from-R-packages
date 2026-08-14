program test_repeat_rng
   use tabu_search, only : dp, i8, tabu_control, tabu_result, run_tabu_search
   use tabu_search_rng, only : tabu_rng
   implicit none

   type(tabu_control) :: control
   type(tabu_result) :: a, b
   type(tabu_rng) :: rng
   integer :: initial(8)
   real(dp) :: u1, u2

   call rng%seed(1_i8)
   u1 = rng%uniform()
   u2 = rng%uniform()
   call check(abs(u1 - 16807.0_dp / 2147483647.0_dp) < 1.0e-15_dp, "Park-Miller first draw")
   call check(abs(u2 - 282475249.0_dp / 2147483647.0_dp) < 1.0e-15_dp, "Park-Miller second draw")

   initial = [1,0,0,1,0,1,0,0]
   control%iters = 9
   control%neigh = 4
   control%list_size = 2
   control%n_restarts = 0
   control%repeat_all = 3
   control%seed = 12345_i8

   call run_tabu_search(8, objective, a, control, initial)
   call run_tabu_search(8, objective, b, control, initial)

   call check(a%n_records == 2 * control%iters * control%repeat_all, "repeat history length")
   call check(all(a%config_keep == b%config_keep), "seeded configuration reproducibility")
   call check(all(abs(a%utility_keep - b%utility_keep) < 1.0e-15_dp), "seeded utility reproducibility")

   print '(a)', "test_repeat_rng: PASS"

contains

   function objective(config) result(value)
      integer, intent(in) :: config(:)
      real(dp) :: value
      integer :: i

      value = 0.0_dp
      do i = 1, size(config)
         value = value + real(i * config(i), dp)
      end do
   end function objective

   subroutine check(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         write(*,'(a)') "FAIL: " // trim(message)
         error stop 1
      end if
   end subroutine check

end program test_repeat_rng
