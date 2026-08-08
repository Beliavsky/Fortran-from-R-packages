! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of computational code from the R package dfoptim.
! Original authors: Ravi Varadhan, Hans W. Borchers, and Vincent Bechard.

module dfoptim_kinds
   implicit none
   private
   integer, parameter, public :: dp = kind(1.0d0)
end module dfoptim_kinds
