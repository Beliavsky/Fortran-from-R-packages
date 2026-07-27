! SPDX-License-Identifier: MIT
! Copyright (c) 2014-2022 Robert L. McDonald
module derivmkts_bonds
    use derivmkts_kinds, only: dp
    implicit none
    private
    public :: bondpv, bondyield, duration, convexity
contains

    pure real(dp) function bondpv(coupon, mat, yield_rate, principal, freq) result(pv)
        real(dp), intent(in) :: coupon, mat, yield_rate
        real(dp), intent(in), optional :: principal
        integer, intent(in), optional :: freq
        real(dp) :: par, coupon_period, yield_period
        integer :: frequency, nperiod, i

        par = 1000.0_dp
        if (present(principal)) par = principal
        frequency = 1
        if (present(freq)) frequency = freq
        nperiod = nint(mat*real(frequency, dp))
        coupon_period = coupon/real(frequency, dp)
        yield_period = yield_rate/real(frequency, dp)
        pv = 0.0_dp
        do i = 1, nperiod
            pv = pv + coupon_period/(1.0_dp + yield_period)**i
        end do
        pv = pv + par/(1.0_dp + yield_period)**nperiod
    end function bondpv

    real(dp) function bondyield(price, coupon, mat, principal, freq, ok) result(yield_rate)
        real(dp), intent(in) :: price, coupon, mat
        real(dp), intent(in), optional :: principal
        integer, intent(in), optional :: freq
        logical, intent(out), optional :: ok
        real(dp) :: lo, hi, mid, flo, fm, par
        integer :: frequency, iter

        par = 1000.0_dp
        if (present(principal)) par = principal
        frequency = 1
        if (present(freq)) frequency = freq
        lo = -0.25_dp
        hi = 20.0_dp
        flo = bondpv(coupon, mat, lo, par, frequency) - price
        do iter = 1, 200
            mid = 0.5_dp*(lo + hi)
            fm = bondpv(coupon, mat, mid, par, frequency) - price
            if (abs(fm) < 1.0e-11_dp .or. hi-lo < 1.0e-12_dp) exit
            if (flo*fm <= 0.0_dp) then
                hi = mid
            else
                lo = mid
                flo = fm
            end if
        end do
        yield_rate = mid
        if (present(ok)) ok = .true.
    end function bondyield

    real(dp) function duration(price, coupon, mat, principal, freq, modified) result(dur)
        real(dp), intent(in) :: price, coupon, mat
        real(dp), intent(in), optional :: principal
        integer, intent(in), optional :: freq
        logical, intent(in), optional :: modified
        real(dp) :: par, yield_rate, yield_period, coupon_period, cashflow
        integer :: frequency, nperiod, i
        logical :: use_modified

        par = 1000.0_dp
        if (present(principal)) par = principal
        frequency = 1
        if (present(freq)) frequency = freq
        use_modified = .false.
        if (present(modified)) use_modified = modified
        yield_rate = bondyield(price, coupon, mat, par, frequency)
        yield_period = yield_rate/real(frequency, dp)
        coupon_period = coupon/real(frequency, dp)
        nperiod = nint(mat*real(frequency, dp))
        dur = 0.0_dp
        do i = 1, nperiod
            cashflow = coupon_period
            if (i == nperiod) cashflow = cashflow + par
            dur = dur + real(i, dp)*cashflow/(1.0_dp + yield_period)**i
        end do
        dur = dur/price/real(frequency, dp)
        if (use_modified) dur = dur/(1.0_dp + yield_period)
    end function duration

    real(dp) function convexity(price, coupon, mat, principal, freq) result(conv)
        real(dp), intent(in) :: price, coupon, mat
        real(dp), intent(in), optional :: principal
        integer, intent(in), optional :: freq
        real(dp) :: par, yield_rate, yield_period, coupon_period, cashflow, weighted_pv
        integer :: frequency, nperiod, i

        par = 1000.0_dp
        if (present(principal)) par = principal
        frequency = 1
        if (present(freq)) frequency = freq
        yield_rate = bondyield(price, coupon, mat, par, frequency)
        yield_period = yield_rate/real(frequency, dp)
        coupon_period = coupon/real(frequency, dp)
        nperiod = nint(mat*real(frequency, dp))
        weighted_pv = 0.0_dp
        do i = 1, nperiod
            cashflow = coupon_period
            if (i == nperiod) cashflow = cashflow + par
            weighted_pv = weighted_pv + (real(i,dp) + real(i*i,dp))*cashflow/ &
                (1.0_dp + yield_period)**i
        end do
        conv = weighted_pv/price/(1.0_dp + yield_period)**2/real(frequency*frequency, dp)
    end function convexity

end module derivmkts_bonds
