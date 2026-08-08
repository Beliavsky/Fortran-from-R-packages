program test_benchmarks
   use pso, only : dp, pso_test_problem, make_test_problem, parabola, &
      griewank, rastrigin, ackley
   implicit none
   type(pso_test_problem) :: prob
   real(dp) :: x2(2), x3(3)

   x3 = 0.0_dp
   if (abs(parabola(x3)) > 1.0e-14_dp) error stop "parabola value mismatch"
   x2 = 0.0_dp
   if (abs(griewank(x2)) > 1.0e-14_dp) error stop "griewank value mismatch"
   if (abs(rastrigin(x2)) > 1.0e-14_dp) error stop "rastrigin value mismatch"
   if (abs(ackley(x2)) > 1.0e-12_dp) error stop "ackley value mismatch"

   call make_test_problem("rastrigin", prob, ntest=3)
   if (prob%n /= 2 .or. prob%maxf /= 3000) error stop "test problem metadata mismatch"
   if (.not. associated(prob%f) .or. .not. associated(prob%grad)) &
      error stop "test problem callbacks not associated"
   print *, "test_benchmarks: PASS"
end program test_benchmarks
