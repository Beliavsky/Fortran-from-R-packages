program test_persistent
   use highs
   implicit none
   type(highs_model) :: model
   type(highs_solver) :: solver, solver2
   type(highs_solution) :: sol
   type(highs_basis) :: basis
   integer(highs_int) :: status
   real(dp) :: a(1,2)
   character(len=*), parameter :: filename = "highs_fortran_test_model.mps"

   if (.not. highs_backend_available()) then
      print *, "test_persistent: SKIP (run scripts/build_backend first)"
      stop
   end if

   a = reshape([1.0_dp,1.0_dp], [1,2])
   call highs_model_from_dense(model, [1.0_dp,1.0_dp], [0.0_dp,0.0_dp], &
      [highs_default_infinity,highs_default_infinity], status, a=a, lhs=[1.0_dp], &
      rhs=[highs_default_infinity])
   call highs_new_solver(solver, status)
   call highs_pass_model(solver, model, status)
   call highs_run(solver, status)
   call highs_get_solution(solver, sol)
   if (abs(sol%objective_value - 1.0_dp) > 1.0e-8_dp) error stop "initial persistent solve failed"

   call highs_change_costs(solver, [1,2], [2.0_dp,1.0_dp], status)
   call highs_run(solver, status)
   call highs_get_solution(solver, sol)
   if (abs(sol%col_value(1)) > 1.0e-7_dp .or. abs(sol%col_value(2)-1.0_dp) > 1.0e-7_dp) &
      error stop "modified objective failed"

   call highs_get_basis(solver, basis, status)
   if (.not. basis%valid) error stop "basis unavailable"
   call highs_set_basis(solver, basis, status)
   if (status == highs_status_error) error stop "basis reset failed"

   call highs_write_model(solver, filename, status)
   if (status == highs_status_error) error stop "model write failed"
   call highs_new_solver(solver2, status)
   call highs_read_model(solver2, filename, status)
   call highs_run(solver2, status)
   call highs_get_solution(solver2, sol)
   if (sol%model_status /= highs_model_optimal .or. abs(sol%objective_value-1.0_dp) > 1.0e-8_dp) &
      error stop "model read failed"

   call highs_destroy_solver(solver)
   call highs_destroy_solver(solver2)
   open(unit=10, file=filename, status="old")
   close(10, status="delete")
   print *, "test_persistent: PASS"
end program test_persistent
