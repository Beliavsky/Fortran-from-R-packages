! SPDX-License-Identifier: Artistic-2.0
module ldhmm_status
   implicit none
   private
   integer, parameter, public :: LDHMM_SUCCESS = 0
   integer, parameter, public :: LDHMM_INVALID_ARGUMENT = 1
   integer, parameter, public :: LDHMM_ALLOCATION_ERROR = 2
   integer, parameter, public :: LDHMM_NUMERICAL_ERROR = 3
   integer, parameter, public :: LDHMM_MAX_ITERATIONS = 4
   integer, parameter, public :: LDHMM_LINE_SEARCH_FAILED = 5
end module ldhmm_status
