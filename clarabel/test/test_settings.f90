program test_settings
   use clarabel
   implicit none
   type(clarabel_settings) :: settings
   type(clarabel_cone) :: cone
   logical :: ok
   character(len=:), allocatable :: message

   settings = default_clarabel_settings()
   if (settings%max_iter /= 200) error stop "max_iter default"
   if (settings%direct_solve_method /= direct_solver_qdldl) error stop "direct solver default"
   if (.not. settings%presolve_enable) error stop "presolve default"
   if (settings%chordal_decomposition_enable) error stop "R wrapper chordal default"
   if (settings%chordal_decomposition_compact) error stop "R wrapper compact default"
   if (settings%chordal_decomposition_complete_dual) error stop "R wrapper complete-dual default"
   if (settings%chordal_decomposition_merge_method /= chordal_merge_none) error stop "merge default"

   cone = generalized_power_cone([0.25_dp, 0.75_dp], 2)
   call cone%validate(ok, message)
   if (.not. ok) error stop message
   cone = generalized_power_cone([1.0_dp], 2)
   call cone%validate(ok, message)
   if (ok) error stop "one-exponent generalized power cone accepted"
   print *, "test_settings: PASS"
end program test_settings
