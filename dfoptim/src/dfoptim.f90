! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of computational code from the R package dfoptim.
! Original authors: Ravi Varadhan, Hans W. Borchers, and Vincent Bechard.

module dfoptim
   use dfoptim_kinds, only : dp
   use dfoptim_interfaces
   use dfoptim_hooke_jeeves, only : hjk, hjkb
   use dfoptim_nelder_mead, only : nmk, nmkb
   use dfoptim_mads, only : mads
   implicit none
   public
end module dfoptim
