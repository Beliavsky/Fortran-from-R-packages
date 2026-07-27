! SPDX-License-Identifier: GPL-3.0-only
! Derived from LSMonteCarlo 1.0 by Mikhail A. Beketov.
! Copyright (C) 2013 Mikhail A. Beketov.
module lsmc_european
    use lsmc_kinds, only : dp
    use lsmc_math, only : normal_cdf
    implicit none
    private

    public :: eu_call_bs
    public :: eu_put_bs

contains

    pure function eu_call_bs(spot, sigma, strike, rate, dividend, maturity) result(value)
        real(dp), intent(in) :: spot
        real(dp), intent(in) :: sigma
        real(dp), intent(in) :: strike
        real(dp), intent(in) :: rate
        real(dp), intent(in) :: dividend
        real(dp), intent(in) :: maturity
        real(dp) :: d1
        real(dp) :: d2
        real(dp) :: value

        call validate_inputs(spot, sigma, strike, maturity)
        if (maturity <= tiny(1.0_dp)) then
            value = max(spot - strike, 0.0_dp)
        else if (sigma <= tiny(1.0_dp)) then
            value = exp(-rate * maturity) * max(spot * exp((rate - dividend) * maturity) - strike, 0.0_dp)
        else
            d1 = (log(spot / strike) + (rate - dividend + 0.5_dp * sigma * sigma) * maturity) / &
                (sigma * sqrt(maturity))
            d2 = d1 - sigma * sqrt(maturity)
            value = spot * exp(-dividend * maturity) * normal_cdf(d1) - &
                strike * exp(-rate * maturity) * normal_cdf(d2)
        end if
    end function eu_call_bs

    pure function eu_put_bs(spot, sigma, strike, rate, dividend, maturity) result(value)
        real(dp), intent(in) :: spot
        real(dp), intent(in) :: sigma
        real(dp), intent(in) :: strike
        real(dp), intent(in) :: rate
        real(dp), intent(in) :: dividend
        real(dp), intent(in) :: maturity
        real(dp) :: d1
        real(dp) :: d2
        real(dp) :: value

        call validate_inputs(spot, sigma, strike, maturity)
        if (maturity <= tiny(1.0_dp)) then
            value = max(strike - spot, 0.0_dp)
        else if (sigma <= tiny(1.0_dp)) then
            value = exp(-rate * maturity) * max(strike - spot * exp((rate - dividend) * maturity), 0.0_dp)
        else
            d1 = (log(spot / strike) + (rate - dividend + 0.5_dp * sigma * sigma) * maturity) / &
                (sigma * sqrt(maturity))
            d2 = d1 - sigma * sqrt(maturity)
            value = strike * exp(-rate * maturity) * normal_cdf(-d2) - &
                spot * exp(-dividend * maturity) * normal_cdf(-d1)
        end if
    end function eu_put_bs

    pure subroutine validate_inputs(spot, sigma, strike, maturity)
        real(dp), intent(in) :: spot
        real(dp), intent(in) :: sigma
        real(dp), intent(in) :: strike
        real(dp), intent(in) :: maturity

        if (spot <= 0.0_dp) error stop "Black-Scholes: spot must be positive"
        if (strike <= 0.0_dp) error stop "Black-Scholes: strike must be positive"
        if (sigma < 0.0_dp) error stop "Black-Scholes: sigma must be nonnegative"
        if (maturity < 0.0_dp) error stop "Black-Scholes: maturity must be nonnegative"
    end subroutine validate_inputs

end module lsmc_european
