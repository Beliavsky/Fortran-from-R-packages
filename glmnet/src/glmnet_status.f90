! SPDX-License-Identifier: GPL-2.0-only
module glmnet_status
   implicit none
   private
   integer, parameter, public :: glmnet_success = 0
   integer, parameter, public :: glmnet_invalid_argument = 1
   integer, parameter, public :: glmnet_nonfinite_input = 2
   integer, parameter, public :: glmnet_all_predictors_constant = 3
   integer, parameter, public :: glmnet_max_iterations = 4
   integer, parameter, public :: glmnet_numerical_failure = 5
   integer, parameter, public :: glmnet_unsupported = 6
   integer, parameter, public :: glmnet_no_events = 7
   integer, parameter, public :: glmnet_invalid_response = 8
   public :: glmnet_status_message
contains
   pure function glmnet_status_message(status) result(message)
      integer, intent(in) :: status
      character(len=:), allocatable :: message
      select case (status)
      case (glmnet_success)
         message = 'success'
      case (glmnet_invalid_argument)
         message = 'invalid argument'
      case (glmnet_nonfinite_input)
         message = 'nonfinite input'
      case (glmnet_all_predictors_constant)
         message = 'all predictors are constant or excluded'
      case (glmnet_max_iterations)
         message = 'iteration limit reached'
      case (glmnet_numerical_failure)
         message = 'numerical failure'
      case (glmnet_unsupported)
         message = 'unsupported option'
      case (glmnet_no_events)
         message = 'no events in survival response'
      case (glmnet_invalid_response)
         message = 'invalid response'
      case default
         message = 'unknown status'
      end select
   end function glmnet_status_message
end module glmnet_status
