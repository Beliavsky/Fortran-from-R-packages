! SPDX-License-Identifier: GPL-2.0-or-later
module maxlik_status
  implicit none
  private

  integer, parameter, public :: MAXLIK_SUCCESS_GRADIENT = 1
  integer, parameter, public :: MAXLIK_SUCCESS_VALUE = 2
  integer, parameter, public :: MAXLIK_STEP_FAILURE = 3
  integer, parameter, public :: MAXLIK_ITERATION_LIMIT = 4
  integer, parameter, public :: MAXLIK_CONSTRAINT_FAILURE = 5
  integer, parameter, public :: MAXLIK_INVALID_START = 100
  integer, parameter, public :: MAXLIK_INVALID_INPUT = 101
  integer, parameter, public :: MAXLIK_EVALUATION_ERROR = 102
  integer, parameter, public :: MAXLIK_SINGULAR_HESSIAN = 103

  public :: maxlik_message

contains

  pure function maxlik_message(code) result(message)
    integer, intent(in) :: code
    character(len=96) :: message

    select case (code)
    case (MAXLIK_SUCCESS_GRADIENT)
      message = 'gradient close to zero'
    case (MAXLIK_SUCCESS_VALUE)
      message = 'successive function values within tolerance'
    case (MAXLIK_STEP_FAILURE)
      message = 'could not find a higher point'
    case (MAXLIK_ITERATION_LIMIT)
      message = 'iteration limit exceeded'
    case (MAXLIK_CONSTRAINT_FAILURE)
      message = 'constraint tolerance was not reached'
    case (MAXLIK_INVALID_START)
      message = 'initial value is invalid or non-finite'
    case (MAXLIK_INVALID_INPUT)
      message = 'invalid problem or control input'
    case (MAXLIK_EVALUATION_ERROR)
      message = 'objective or derivative evaluation failed'
    case (MAXLIK_SINGULAR_HESSIAN)
      message = 'Hessian or information matrix is singular'
    case default
      write(message, '(a,i0)') 'unknown return code ', code
    end select
  end function maxlik_message

end module maxlik_status
