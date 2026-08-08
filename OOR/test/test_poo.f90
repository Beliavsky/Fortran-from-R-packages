module test_poo_callbacks
   use oor, only : dp
   implicit none
contains
   function objective(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      y = -(x - 0.35_dp)**2
   end function objective
end module test_poo_callbacks

program test_poo
   use oor, only : dp, poo, poo_result, set_random_seed
   use test_poo_callbacks, only : objective
   implicit none
   type(poo_result) :: result

   call set_random_seed(12)
   call poo(objective, 300, 0.0_dp, result, rhomax=10, nu=1.0_dp)
   if (abs(result%par - 0.35_dp) > 5.0e-4_dp) error stop "POO parameter failure"
   if (result%value < -1.0e-6_dp) error stop "POO value failure"
   if (result%evaluations < 300) error stop "POO evaluation accounting failure"
   if (size(result%tree%leaves) < 2) error stop "POO tree failure"
   print *, "test_poo: PASS"
end program test_poo
