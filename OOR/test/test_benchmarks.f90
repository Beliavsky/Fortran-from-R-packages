module test_benchmark_callbacks
   use oor, only : dp, guirland, sin1, difficult, difficult2, double_sine
   implicit none
contains
   subroutine run_checks()
      real(dp), parameter :: tol = 2.0e-12_dp
      call check(abs(guirland(0.25_dp) - 0.5987992001326592_dp) < tol, "guirland")
      call check(abs(sin1(0.25_dp) - 0.4756537104464140_dp) < tol, "sin1")
      call check(abs(difficult(0.25_dp) - 0.9200056958555479_dp) < tol, "difficult")
      call check(abs(difficult2(0.25_dp) + 0.0625_dp) < tol, "difficult2")
      call check(abs(double_sine(0.25_dp) + 0.55_dp) < tol, "double_sine")
      call check(abs(double_sine(0.5_dp)) < tol, "double_sine center")
   end subroutine run_checks

   subroutine check(ok, label)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: label
      if (.not. ok) error stop "benchmark failure: " // label
   end subroutine check
end module test_benchmark_callbacks

program test_benchmarks
   use test_benchmark_callbacks, only : run_checks
   implicit none
   call run_checks()
   print *, "test_benchmarks: PASS"
end program test_benchmarks
