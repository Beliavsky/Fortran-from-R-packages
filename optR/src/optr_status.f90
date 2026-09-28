module optr_status
   implicit none
   private

   integer, parameter, public :: optr_success = 0
   integer, parameter, public :: optr_invalid_shape = 1
   integer, parameter, public :: optr_singular = 2
   integer, parameter, public :: optr_not_positive_definite = 3
   integer, parameter, public :: optr_no_convergence = 4
   integer, parameter, public :: optr_invalid_argument = 5

end module optr_status
