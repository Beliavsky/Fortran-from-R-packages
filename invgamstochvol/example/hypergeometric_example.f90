! SPDX-License-Identifier: MIT
program hypergeometric_example
   use invgamstochvol
   implicit none

   real(dp) :: value
   integer :: status

   value = ourgeo(1.5_dp, 1.9_dp, 1.2_dp, 0.7_dp, status=status)
   write (*, '(a,f18.10)') '2F1(1.5, 1.9; 1.2; 0.7) = ', value
   write (*, '(a,i0)') 'status = ', status
end program hypergeometric_example
