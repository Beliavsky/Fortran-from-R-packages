module e1071_constants
    use e1071_kinds, only: dp
    implicit none
    private

    real(dp), parameter, public :: e1071_pi = acos(-1.0_dp)
    real(dp), parameter, public :: e1071_tiny = sqrt(tiny(1.0_dp))
end module e1071_constants
