module circstats_kinds
    implicit none
    integer, parameter, public :: dp = kind(1.0d0)
    real(dp), parameter, public :: pi = acos(-1.0_dp)
    real(dp), parameter, public :: twopi = 2.0_dp*pi
end module circstats_kinds
