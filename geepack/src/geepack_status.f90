! Modern Fortran translation of computational code from geepack 1.3-13.
! Upstream license: GPL (>= 3). See LICENSE, NOTICE.md, and PROVENANCE.md.
! Exact upstream copyright-holder attribution is preserved in NOTICE.md and upstream/DESCRIPTION.
module geepack_status
   implicit none
   private

   integer, parameter, public :: GEE_OK = 0
   integer, parameter, public :: GEE_ERR_SHAPE = 1
   integer, parameter, public :: GEE_ERR_ARGUMENT = 2
   integer, parameter, public :: GEE_ERR_SINGULAR = 3
   integer, parameter, public :: GEE_ERR_MAXITER = 4
   integer, parameter, public :: GEE_ERR_INVALID_MEAN = 5
   integer, parameter, public :: GEE_ERR_CORRELATION = 6

end module geepack_status
