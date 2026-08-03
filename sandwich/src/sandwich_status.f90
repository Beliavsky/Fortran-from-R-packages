! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module sandwich_status
   implicit none
   private
   public :: set_status

   integer, parameter, public :: SANDWICH_SUCCESS = 0
   integer, parameter, public :: SANDWICH_INVALID_ARGUMENT = 1
   integer, parameter, public :: SANDWICH_DIMENSION_MISMATCH = 2
   integer, parameter, public :: SANDWICH_SINGULAR_MATRIX = 3
   integer, parameter, public :: SANDWICH_INSUFFICIENT_DATA = 4
   integer, parameter, public :: SANDWICH_NUMERICAL_FAILURE = 5
   integer, parameter, public :: SANDWICH_UNSUPPORTED = 6

contains

   subroutine set_status(status, value)
      integer, intent(out), optional :: status
      integer, intent(in) :: value
      if (present(status)) status = value
   end subroutine set_status

end module sandwich_status
