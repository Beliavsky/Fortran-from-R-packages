! SPDX-License-Identifier: MIT
! SPDX-FileComment: Result containers for the Fortran stringr translation.
module stringr_types
   !! Defines rectangular results needed by split and location operations.
   implicit none
   private

   type, public :: split_result_type
      !! Stores split fields by input row and the valid field count per row.
      character(len=:), allocatable :: values(:,:)
      integer, allocatable :: counts(:)
   end type split_result_type

   type, public :: location_type
      !! Stores one-based inclusive match bounds, with zeros for no match.
      integer :: start = 0
      integer :: end = 0
   end type location_type

end module stringr_types
