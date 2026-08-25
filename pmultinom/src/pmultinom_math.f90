! Numerical helpers for pmultinom-fortran.
module pmultinom_math
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use pmultinom_kinds, only : dp
    implicit none
    private

    public :: log_poisson_pmf, integer_support

contains

    pure real(dp) function log_poisson_pmf(k, lambda) result(ans)
        integer, intent(in) :: k
        real(dp), intent(in) :: lambda

        if (k < 0 .or. lambda < 0.0_dp) then
            ans = -huge(1.0_dp)
        else if (lambda == 0.0_dp) then
            if (k == 0) then
                ans = 0.0_dp
            else
                ans = -huge(1.0_dp)
            end if
        else
            ans = -lambda + real(k, dp) * log(lambda) - log_gamma(real(k + 1, dp))
        end if
    end function log_poisson_pmf

    pure subroutine integer_support(lower, upper, nmax, ilo, ihi)
        real(dp), intent(in) :: lower, upper
        integer, intent(in) :: nmax
        integer, intent(out) :: ilo, ihi
        real(dp) :: fl

        if (.not. ieee_is_finite(lower)) then
            if (lower < 0.0_dp) then
                ilo = 0
            else
                ilo = nmax + 1
            end if
        else if (lower < -real(huge(1) - 2, dp)) then
            ilo = 0
        else if (lower >= real(nmax, dp)) then
            ilo = nmax + 1
        else
            fl = floor(lower)
            ilo = max(0, int(fl) + 1)
        end if

        if (.not. ieee_is_finite(upper)) then
            if (upper > 0.0_dp) then
                ihi = nmax
            else
                ihi = -1
            end if
        else if (upper < 0.0_dp) then
            ihi = -1
        else if (upper >= real(nmax, dp)) then
            ihi = nmax
        else
            ihi = min(nmax, int(floor(upper)))
        end if
    end subroutine integer_support

end module pmultinom_math
