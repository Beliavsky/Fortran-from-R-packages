! SPDX-License-Identifier: MIT
! SPDX-FileComment: Public facade for the Fortran stringr translation.
module stringr
   !! Re-exports the portable fixed-string subset of stringr.
   use stringr_scalar
   use stringr_types
   use stringr_vector
   implicit none
   public
end module stringr
