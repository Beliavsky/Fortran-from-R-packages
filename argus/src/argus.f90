! Public umbrella module for argus-fortran.
! SPDX-License-Identifier: GPL-2.0-or-later
module argus
   use argus_kinds, only : dp
   use argus_distribution, only : dargus, pargus, qargus, &
      dargus_recycle, pargus_recycle, qargus_recycle, &
      rargus, rargus_varying, seed_argus_rng, ARGUS_INVERSION, ARGUS_ROU
   implicit none
   public
end module argus
