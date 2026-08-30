module proxy_ieee
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
    use proxy_kinds, only: dp
    implicit none
    private

    public :: proxy_nan, proxy_is_missing

contains

    pure function proxy_nan() result(value)
        real(dp) :: value

        value = ieee_value(0.0_dp, ieee_quiet_nan)
    end function proxy_nan

    pure elemental function proxy_is_missing(value) result(missing)
        real(dp), intent(in) :: value !! Numeric value; NaN denotes a missing observation, matching proxy's NA/NaN exclusion
        !! semantics.
        logical :: missing

        missing = ieee_is_nan(value)
    end function proxy_is_missing

end module proxy_ieee
