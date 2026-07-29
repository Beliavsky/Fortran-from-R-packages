! SPDX-License-Identifier: GPL-2.0-only
!
! Public umbrella module for the GCPM modern Fortran translation.
module gcpm
   use gcpm_kinds
   use gcpm_types
   use gcpm_math
   use gcpm_portfolio
   use gcpm_analytical
   use gcpm_simulation
   use gcpm_risk
   use gcpm_csv
   implicit none
   public
end module gcpm
