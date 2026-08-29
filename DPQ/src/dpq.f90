! SPDX-License-Identifier: GPL-3.0-or-later
!
! Umbrella module for the DPQ modern Fortran port.
! The combined distribution is GPL-3.0-or-later because it contains the
! GPL-3.0-or-later TOMS 1006 component. Individual source files retain their
! more permissive SPDX identifiers where applicable.
module dpq
   use r_compat, only: dp
   use dpq_core
   use dpq_gamma_discrete
   use dpq_normal_beta
   use dpq_hyper
   use dpq_nchisq
   use dpq_t
   use dpq_wiener
   use dpq_toms1006
   implicit none
   public
end module dpq
