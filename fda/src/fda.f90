! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from the computational code of R package fda 6.3.0.
module fda
   use r_kinds, only : dp
   use fda_basis
   use fda_numeric
   use fda_fd
   use fda_smoothing
   use fda_analysis
   use fda_ode
   implicit none
   public
end module fda
