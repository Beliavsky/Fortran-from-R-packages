program test_legacy_and_odd
   use metaheuristic_opt, only : dp, mh_control, mh_result, ffa, clonalg, goa
   implicit none
   type(mh_control) :: c
   type(mh_result) :: legacy, fixed, odd
   real(dp) :: lo4(4), hi4(4), lo3(3), hi3(3)

   lo4 = -5.0_dp
   hi4 = 5.0_dp
   c%num_population = 24
   c%max_iter = 20
   c%seed = 44
   c%legacy_quirks = .true.
   call ffa(sphere_local, lo4, hi4, legacy, c)
   c%legacy_quirks = .false.
   c%seed = 44
   call ffa(sphere_local, lo4, hi4, fixed, c)
   if (legacy%value < fixed%value - 1.0e-12_dp) error stop 'FFA legacy which.max behavior changed'

   c%legacy_quirks = .true.
   c%selection_size = 6
   c%seed = 55
   call clonalg(sphere_local, lo4, hi4, legacy, c)
   c%legacy_quirks = .false.
   c%seed = 55
   call clonalg(sphere_local, lo4, hi4, fixed, c)
   if (legacy%value < fixed%value - 1.0e-12_dp) error stop 'CLONALG legacy return behavior changed'

   lo3 = -3.0_dp
   hi3 = 3.0_dp
   c%max_iter = 20
   c%seed = 66
   call goa(sphere_local, lo3, hi3, odd, c)
   if (size(odd%par) /= 3) error stop 'GOA odd-dimensional result has wrong size'
   if (any(odd%par < lo3) .or. any(odd%par > hi3)) error stop 'GOA odd-dimensional result out of bounds'
contains
   real(dp) function sphere_local(x) result(f)
      real(dp), intent(in) :: x(:)
      f = sum(x*x)
   end function sphere_local
end program test_legacy_and_odd
