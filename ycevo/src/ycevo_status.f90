module ycevo_status
   implicit none
   private

   integer, parameter, public :: ycevo_success = 0
   integer, parameter, public :: ycevo_err_input = 1
   integer, parameter, public :: ycevo_err_singular = 2
   integer, parameter, public :: ycevo_err_allocation = 3
end module ycevo_status
