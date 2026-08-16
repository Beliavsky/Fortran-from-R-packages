module compoissonreg_kinds
   implicit none
   private
   integer, parameter, public :: dp = kind(1.0d0)
   real(dp), parameter, public :: pi_dp = acos(-1.0_dp)
end module compoissonreg_kinds
