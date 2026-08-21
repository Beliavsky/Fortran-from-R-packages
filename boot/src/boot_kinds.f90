module boot_kinds
    implicit none
    private
    integer, parameter, public :: dp = kind(1.0d0)
    real(dp), parameter, public :: pi = acos(-1.0_dp)
end module boot_kinds
