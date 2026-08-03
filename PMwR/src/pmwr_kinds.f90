module pmwr_kinds
   implicit none
   private
   integer, parameter, public :: dp = kind(1.0d0)
   integer, parameter, public :: i8 = selected_int_kind(18)
end module pmwr_kinds
