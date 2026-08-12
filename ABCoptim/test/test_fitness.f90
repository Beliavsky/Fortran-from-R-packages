program test_fitness
  use abcoptim, only : dp, calculate_fitness
  implicit none

  call check(abs(calculate_fitness(0.0_dp) - 1.0_dp) < 1.0e-15_dp, "fitness at zero")
  call check(abs(calculate_fitness(3.0_dp) - 0.25_dp) < 1.0e-15_dp, "positive fitness")
  call check(abs(calculate_fitness(-3.0_dp) - 4.0_dp) < 1.0e-15_dp, "negative fitness")
  print *, "test_fitness: PASS"

contains

  subroutine check(ok, message)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: message
    if (.not. ok) then
      print *, "FAIL: ", trim(message)
      error stop 1
    end if
  end subroutine check

end program test_fitness
