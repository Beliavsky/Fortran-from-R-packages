! Experimental modern Fortran translation of computational routines from
! the R package tseries 0.10-62. Original authors include Adrian Trapletti
! and Kurt Hornik; Blake LeBaron contributed the original BDS code.
! Licensed under GPL-2.0-only OR GPL-3.0-only. See LICENSE and NOTICE.

module tseries_kinds
   implicit none
   private
   public :: dp
   integer, parameter :: dp = kind(1.0d0)
end module tseries_kinds
