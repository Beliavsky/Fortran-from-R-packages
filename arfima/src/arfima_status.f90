module arfima_status
  implicit none
  private
  integer, parameter, public :: arfima_ok = 0
  integer, parameter, public :: arfima_invalid_input = 1
  integer, parameter, public :: arfima_not_stationary = 2
  integer, parameter, public :: arfima_not_positive_definite = 3
  integer, parameter, public :: arfima_singular = 4
  integer, parameter, public :: arfima_no_convergence = 5
  integer, parameter, public :: arfima_allocation_error = 6
  integer, parameter, public :: arfima_numerical_error = 7
end module arfima_status
