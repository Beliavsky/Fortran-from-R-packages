! SPDX-License-Identifier: GPL-2.0-or-later
!
! Modern Fortran translation of the computational algorithms in the R package
! corpmetrics 1.0 by Pavlos Pantatosakis.
module corpmetrics
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    implicit none
    private

    integer, parameter, public :: dp = kind(1.0d0)

    integer, parameter, public :: cm_success = 0
    integer, parameter, public :: cm_invalid_input = 1
    integer, parameter, public :: cm_dimension_mismatch = 2
    integer, parameter, public :: cm_singular = 3
    integer, parameter, public :: cm_root_not_bracketed = 4
    integer, parameter, public :: cm_root_failure = 5

    integer, parameter, public :: ddm_zero_growth = 1
    integer, parameter, public :: ddm_gordon = 2
    integer, parameter, public :: ddm_differential = 3

    integer, parameter, public :: capm_undervalued = -1
    integer, parameter, public :: capm_fairly_valued = 0
    integer, parameter, public :: capm_overvalued = 1

    type, public :: balance_sheet_result
        real(dp) :: total_assets = 0.0_dp
        real(dp) :: total_liabilities = 0.0_dp
        real(dp) :: equity = 0.0_dp
        real(dp) :: working_capital = 0.0_dp
        real(dp) :: current_ratio = 0.0_dp
        real(dp) :: acid_test_ratio = 0.0_dp
        real(dp) :: leverage_ratio = 0.0_dp
        real(dp) :: debt_to_equity_ratio = 0.0_dp
    end type balance_sheet_result

    type, public :: capm_result
        real(dp) :: expected_asset_return = 0.0_dp
        real(dp) :: expected_market_return = 0.0_dp
        real(dp) :: risk_free_rate = 0.0_dp
        real(dp) :: beta = 0.0_dp
        real(dp) :: required_return = 0.0_dp
        integer :: valuation = capm_fairly_valued
    end type capm_result

    type, public :: ddm_result
        integer :: model = 0
        real(dp) :: value = 0.0_dp
    end type ddm_result

    type, public :: fixed_income_result
        real(dp) :: price = 0.0_dp
        real(dp) :: macaulay_duration = 0.0_dp
        real(dp) :: modified_duration = 0.0_dp
        integer :: periods = 0
        logical :: semiannual = .false.
    end type fixed_income_result

    type, public :: investment_result
        real(dp) :: npv = 0.0_dp
        real(dp) :: irr = 0.0_dp
        real(dp) :: irr_percent = 0.0_dp
    end type investment_result

    type, public :: income_statement_result
        real(dp) :: gross_profit_margin_percent = 0.0_dp
        real(dp) :: net_profit_margin_percent = 0.0_dp
        real(dp) :: earnings_per_share = 0.0_dp
        real(dp) :: price_earnings_ratio = 0.0_dp
        logical :: has_eps = .false.
        logical :: has_pe = .false.
    end type income_statement_result

    type, public :: loan_result
        real(dp) :: installment = 0.0_dp
        real(dp) :: total_repayment = 0.0_dp
        integer, allocatable :: period(:)
        real(dp), allocatable :: interest(:)
        real(dp), allocatable :: principal(:)
        real(dp), allocatable :: balance(:)
    end type loan_result

    public :: balsh, capm, ddm, fis, idm, insta, loan
    public :: net_present_value, internal_rate_of_return
    public :: valuation_label, ddm_model_label, round2_r

contains

    subroutine set_status(status, value)
        integer, intent(out), optional :: status
        integer, intent(in) :: value

        if (present(status)) status = value
    end subroutine set_status

    pure logical function near_zero(x, scale) result(is_zero)
        real(dp), intent(in) :: x, scale

        is_zero = abs(x) <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(scale))
    end function near_zero

    pure real(dp) function nan_dp() result(value)
        value = ieee_value(0.0_dp, ieee_quiet_nan)
    end function nan_dp

    pure real(dp) function mean_dp(x) result(value)
        real(dp), intent(in) :: x(:)

        value = sum(x) / real(size(x), dp)
    end function mean_dp

    pure real(dp) function sample_variance(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: xbar

        xbar = mean_dp(x)
        value = sum((x - xbar) ** 2) / real(size(x) - 1, dp)
    end function sample_variance

    pure real(dp) function sample_covariance(x, y) result(value)
        real(dp), intent(in) :: x(:), y(:)
        real(dp) :: xbar, ybar

        xbar = mean_dp(x)
        ybar = mean_dp(y)
        value = sum((x - xbar) * (y - ybar)) / real(size(x) - 1, dp)
    end function sample_covariance

    pure real(dp) function round2_r(x) result(value)
        real(dp), intent(in) :: x
        real(dp) :: scaled, lower, fraction
        integer :: ilower

        scaled = x * 100.0_dp
        lower = floor(scaled)
        ilower = int(lower)
        fraction = scaled - lower

        if (fraction < 0.5_dp) then
            value = lower / 100.0_dp
        else if (fraction > 0.5_dp) then
            value = (lower + 1.0_dp) / 100.0_dp
        else
            if (mod(abs(ilower), 2) == 0) then
                value = lower / 100.0_dp
            else
                value = (lower + 1.0_dp) / 100.0_dp
            end if
        end if
    end function round2_r

    pure function valuation_label(code) result(label)
        integer, intent(in) :: code
        character(len=:), allocatable :: label

        select case (code)
        case (capm_undervalued)
            label = 'undervalued'
        case (capm_overvalued)
            label = 'overvalued'
        case default
            label = 'fairly valued'
        end select
    end function valuation_label

    pure function ddm_model_label(code) result(label)
        integer, intent(in) :: code
        character(len=:), allocatable :: label

        select case (code)
        case (ddm_zero_growth)
            label = 'Zero Growth Model'
        case (ddm_gordon)
            label = "Gordon's Model"
        case (ddm_differential)
            label = 'Differential Growth Model'
        case default
            label = 'Unknown Model'
        end select
    end function ddm_model_label

    subroutine balsh(fa, ca, inv, fl, cl, result, status)
        real(dp), intent(in) :: fa, ca, inv, fl, cl
        type(balance_sheet_result), intent(out) :: result
        integer, intent(out), optional :: status

        result = balance_sheet_result()
        if (.not. all(ieee_is_finite([fa, ca, inv, fl, cl]))) then
            call set_status(status, cm_invalid_input)
            return
        end if

        result%total_assets = fa + ca
        result%total_liabilities = fl + cl
        result%equity = result%total_assets - result%total_liabilities
        result%working_capital = ca - cl

        if (near_zero(cl, ca) .or. near_zero(result%total_assets, fa + ca) .or. &
            near_zero(result%equity, result%total_assets)) then
            result%current_ratio = nan_dp()
            result%acid_test_ratio = nan_dp()
            result%leverage_ratio = nan_dp()
            result%debt_to_equity_ratio = nan_dp()
            call set_status(status, cm_singular)
            return
        end if

        result%current_ratio = ca / cl
        result%acid_test_ratio = (ca - inv) / cl
        result%leverage_ratio = result%total_liabilities / result%total_assets
        result%debt_to_equity_ratio = result%total_liabilities / result%equity
        call set_status(status, cm_success)
    end subroutine balsh

    subroutine capm(rf, ri, rm, result, status)
        real(dp), intent(in) :: rf(:), ri(:), rm(:)
        type(capm_result), intent(out) :: result
        integer, intent(out), optional :: status
        real(dp) :: market_variance

        result = capm_result()
        if (size(rf) == 0 .or. size(ri) < 2 .or. size(rm) < 2) then
            call set_status(status, cm_invalid_input)
            return
        end if
        if (size(ri) /= size(rm)) then
            call set_status(status, cm_dimension_mismatch)
            return
        end if
        if (.not. all(ieee_is_finite(rf)) .or. .not. all(ieee_is_finite(ri)) .or. &
            .not. all(ieee_is_finite(rm))) then
            call set_status(status, cm_invalid_input)
            return
        end if

        result%expected_asset_return = mean_dp(ri)
        result%expected_market_return = mean_dp(rm)
        result%risk_free_rate = mean_dp(rf)
        market_variance = sample_variance(rm)
        if (near_zero(market_variance, maxval(abs(rm)) ** 2)) then
            result%beta = nan_dp()
            result%required_return = nan_dp()
            call set_status(status, cm_singular)
            return
        end if

        result%beta = sample_covariance(ri, rm) / market_variance
        result%required_return = result%risk_free_rate + result%beta * &
            (result%expected_market_return - result%risk_free_rate)

        if (result%required_return < result%expected_asset_return) then
            result%valuation = capm_undervalued
        else if (result%required_return > result%expected_asset_return) then
            result%valuation = capm_overvalued
        else
            result%valuation = capm_fairly_valued
        end if
        call set_status(status, cm_success)
    end subroutine capm

    subroutine ddm(dividend, required_return, result, g1, g2, period, status)
        real(dp), intent(in) :: dividend, required_return
        type(ddm_result), intent(out) :: result
        real(dp), intent(in), optional :: g1, g2
        integer, intent(in), optional :: period
        integer, intent(out), optional :: status
        real(dp) :: div1, value1, div_n1, value2

        result = ddm_result()
        if (.not. ieee_is_finite(dividend) .or. .not. ieee_is_finite(required_return)) then
            call set_status(status, cm_invalid_input)
            return
        end if

        if (.not. present(g1) .and. .not. present(g2) .and. .not. present(period)) then
            if (near_zero(required_return, dividend)) then
                result%value = nan_dp()
                call set_status(status, cm_singular)
                return
            end if
            result%model = ddm_zero_growth
            result%value = dividend / required_return
        else if (present(g1) .and. .not. present(g2) .and. .not. present(period)) then
            if (.not. ieee_is_finite(g1) .or. near_zero(required_return - g1, required_return)) then
                result%value = nan_dp()
                call set_status(status, cm_singular)
                return
            end if
            result%model = ddm_gordon
            div1 = dividend * (1.0_dp + g1)
            result%value = div1 / (required_return - g1)
        else if (present(g1) .and. present(g2) .and. present(period)) then
            if (.not. ieee_is_finite(g1) .or. .not. ieee_is_finite(g2) .or. period < 0 .or. &
                near_zero(required_return - g1, required_return) .or. &
                near_zero(required_return - g2, required_return) .or. &
                near_zero(1.0_dp + required_return, required_return)) then
                result%value = nan_dp()
                call set_status(status, cm_invalid_input)
                return
            end if
            result%model = ddm_differential
            div1 = dividend * (1.0_dp + g1)
            value1 = div1 / (required_return - g1) * &
                (1.0_dp - (1.0_dp + g1) ** period / (1.0_dp + required_return) ** period)
            div_n1 = dividend * (1.0_dp + g1) ** (period + 1)
            value2 = div_n1 / (required_return - g2) / (1.0_dp + required_return) ** period
            result%value = value1 + value2
        else
            call set_status(status, cm_invalid_input)
            return
        end if

        if (.not. ieee_is_finite(result%value)) then
            call set_status(status, cm_singular)
        else
            call set_status(status, cm_success)
        end if
    end subroutine ddm

    subroutine fis(face_value, coupon_rate, ytm, maturity, result, semiannual, status)
        real(dp), intent(in) :: face_value, coupon_rate, ytm
        integer, intent(in) :: maturity
        type(fixed_income_result), intent(out) :: result
        logical, intent(in), optional :: semiannual
        integer, intent(out), optional :: status
        integer :: t, n
        real(dp) :: periodic_coupon_rate, periodic_yield, coupon, discount, weighted
        logical :: use_semiannual

        result = fixed_income_result()
        use_semiannual = .false.
        if (present(semiannual)) use_semiannual = semiannual

        if (.not. all(ieee_is_finite([face_value, coupon_rate, ytm])) .or. maturity <= 0) then
            call set_status(status, cm_invalid_input)
            return
        end if

        if (use_semiannual) then
            n = 2 * maturity
            periodic_coupon_rate = coupon_rate / 2.0_dp
            periodic_yield = ytm / 2.0_dp
        else
            n = maturity
            periodic_coupon_rate = coupon_rate
            periodic_yield = ytm
        end if
        if (1.0_dp + periodic_yield <= 0.0_dp .or. near_zero(1.0_dp + ytm, ytm)) then
            call set_status(status, cm_invalid_input)
            return
        end if

        coupon = face_value * periodic_coupon_rate
        result%price = 0.0_dp
        weighted = 0.0_dp
        do t = 1, n
            discount = (1.0_dp + periodic_yield) ** t
            result%price = result%price + coupon / discount
            weighted = weighted + real(t, dp) * coupon / discount
        end do
        discount = (1.0_dp + periodic_yield) ** n
        result%price = result%price + face_value / discount
        weighted = weighted + real(n, dp) * face_value / discount

        if (near_zero(result%price, face_value)) then
            call set_status(status, cm_singular)
            return
        end if
        result%macaulay_duration = weighted / result%price
        if (use_semiannual) result%macaulay_duration = result%macaulay_duration / 2.0_dp

        ! This intentionally follows the R source. For semiannual bonds, many
        ! textbooks divide by (1 + ytm/2), while corpmetrics uses (1 + ytm).
        result%modified_duration = result%macaulay_duration / (1.0_dp + ytm)
        result%periods = n
        result%semiannual = use_semiannual
        call set_status(status, cm_success)
    end subroutine fis

    pure real(dp) function net_present_value(cash_flows, costs) result(value)
        real(dp), intent(in) :: cash_flows(:), costs(:)
        integer :: i

        if (size(cash_flows) == 0 .or. size(cash_flows) /= size(costs)) then
            value = nan_dp()
            return
        end if
        value = 0.0_dp
        do i = 1, size(cash_flows)
            if (1.0_dp + costs(i) <= 0.0_dp .and. i > 1) then
                value = nan_dp()
                return
            end if
            value = value + cash_flows(i) / (1.0_dp + costs(i)) ** (i - 1)
        end do
    end function net_present_value

    pure real(dp) function npv_constant_rate(cash_flows, rate) result(value)
        real(dp), intent(in) :: cash_flows(:), rate
        integer :: i

        if (1.0_dp + rate <= 0.0_dp) then
            value = nan_dp()
            return
        end if
        value = 0.0_dp
        do i = 1, size(cash_flows)
            value = value + cash_flows(i) / (1.0_dp + rate) ** (i - 1)
        end do
    end function npv_constant_rate

    subroutine internal_rate_of_return(cash_flows, irr, status, lower, upper, tolerance, max_iterations)
        real(dp), intent(in) :: cash_flows(:)
        real(dp), intent(out) :: irr
        integer, intent(out), optional :: status
        real(dp), intent(in), optional :: lower, upper, tolerance
        integer, intent(in), optional :: max_iterations
        real(dp) :: a, b, c, fa, fb, fc, tol
        integer :: iter, max_iter

        irr = nan_dp()
        if (size(cash_flows) < 2 .or. .not. all(ieee_is_finite(cash_flows))) then
            call set_status(status, cm_invalid_input)
            return
        end if

        a = -1.0_dp + 1.0e-12_dp
        b = 1.0_dp
        if (present(lower)) a = lower
        if (present(upper)) b = upper
        tol = 1.0e-10_dp
        if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
        max_iter = 200
        if (present(max_iterations)) max_iter = max_iterations

        if (a <= -1.0_dp .or. b <= a .or. max_iter <= 0) then
            call set_status(status, cm_invalid_input)
            return
        end if
        fa = npv_constant_rate(cash_flows, a)
        fb = npv_constant_rate(cash_flows, b)
        if (.not. ieee_is_finite(fa) .or. .not. ieee_is_finite(fb)) then
            call set_status(status, cm_root_failure)
            return
        end if
        if (abs(fa) <= tol) then
            irr = a
            call set_status(status, cm_success)
            return
        end if
        if (abs(fb) <= tol) then
            irr = b
            call set_status(status, cm_success)
            return
        end if
        if (fa * fb > 0.0_dp) then
            call set_status(status, cm_root_not_bracketed)
            return
        end if

        do iter = 1, max_iter
            c = 0.5_dp * (a + b)
            fc = npv_constant_rate(cash_flows, c)
            if (.not. ieee_is_finite(fc)) then
                call set_status(status, cm_root_failure)
                return
            end if
            if (abs(fc) <= tol .or. 0.5_dp * abs(b - a) <= tol * max(1.0_dp, abs(c))) then
                irr = c
                call set_status(status, cm_success)
                return
            end if
            if (fa * fc < 0.0_dp) then
                b = c
                fb = fc
            else
                a = c
                fa = fc
            end if
        end do

        irr = 0.5_dp * (a + b)
        call set_status(status, cm_root_failure)
    end subroutine internal_rate_of_return

    subroutine idm(cash_flows, costs, result, status)
        real(dp), intent(in) :: cash_flows(:), costs(:)
        type(investment_result), intent(out) :: result
        integer, intent(out), optional :: status
        integer :: root_status

        result = investment_result()
        if (size(cash_flows) == 0 .or. size(cash_flows) /= size(costs)) then
            call set_status(status, cm_dimension_mismatch)
            return
        end if
        if (.not. all(ieee_is_finite(cash_flows)) .or. .not. all(ieee_is_finite(costs))) then
            call set_status(status, cm_invalid_input)
            return
        end if

        result%npv = net_present_value(cash_flows, costs)
        if (.not. ieee_is_finite(result%npv)) then
            call set_status(status, cm_invalid_input)
            return
        end if
        call internal_rate_of_return(cash_flows, result%irr, root_status)
        if (root_status /= cm_success) then
            result%irr_percent = nan_dp()
            call set_status(status, root_status)
            return
        end if
        result%irr_percent = 100.0_dp * result%irr
        call set_status(status, cm_success)
    end subroutine idm

    subroutine insta(revenue, cost_of_sales, net_income, result, preferred_dividend, shares, &
        price_per_share, status)
        real(dp), intent(in) :: revenue, cost_of_sales, net_income
        type(income_statement_result), intent(out) :: result
        real(dp), intent(in), optional :: preferred_dividend, shares, price_per_share
        integer, intent(out), optional :: status

        result = income_statement_result()
        if (.not. all(ieee_is_finite([revenue, cost_of_sales, net_income])) .or. near_zero(revenue, cost_of_sales)) then
            call set_status(status, cm_invalid_input)
            return
        end if

        result%gross_profit_margin_percent = ((revenue - cost_of_sales) / revenue) * 100.0_dp
        result%net_profit_margin_percent = (net_income / revenue) * 100.0_dp

        if (present(preferred_dividend) .and. present(shares)) then
            if (.not. ieee_is_finite(preferred_dividend) .or. .not. ieee_is_finite(shares) .or. near_zero(shares, net_income)) then
                call set_status(status, cm_invalid_input)
                return
            end if
            result%earnings_per_share = (net_income - preferred_dividend) / shares
            result%has_eps = .true.
        end if

        if (result%has_eps .and. present(price_per_share)) then
            if (.not. ieee_is_finite(price_per_share) .or. near_zero(result%earnings_per_share, net_income)) then
                call set_status(status, cm_singular)
                return
            end if
            result%price_earnings_ratio = price_per_share / result%earnings_per_share
            result%has_pe = .true.
        end if
        call set_status(status, cm_success)
    end subroutine insta

    subroutine loan(amount, rate, periods, result, status)
        real(dp), intent(in) :: amount, rate
        integer, intent(in) :: periods
        type(loan_result), intent(out) :: result
        integer, intent(out), optional :: status
        real(dp) :: loan_balance, interest_payment, principal_payment
        integer :: i

        result%installment = 0.0_dp
        result%total_repayment = 0.0_dp
        if (allocated(result%period)) deallocate(result%period)
        if (allocated(result%interest)) deallocate(result%interest)
        if (allocated(result%principal)) deallocate(result%principal)
        if (allocated(result%balance)) deallocate(result%balance)

        if (.not. ieee_is_finite(amount) .or. .not. ieee_is_finite(rate) .or. periods <= 0 .or. &
            amount < 0.0_dp .or. 1.0_dp + rate <= 0.0_dp) then
            call set_status(status, cm_invalid_input)
            return
        end if

        if (near_zero(rate, 1.0_dp)) then
            result%installment = amount / real(periods, dp)
        else
            result%installment = amount * rate / (1.0_dp - (1.0_dp + rate) ** (-periods))
        end if
        result%total_repayment = result%installment * real(periods, dp)

        allocate(result%period(periods), result%interest(periods), result%principal(periods), &
            result%balance(periods))
        loan_balance = amount
        do i = 1, periods
            interest_payment = loan_balance * rate
            principal_payment = result%installment - interest_payment
            loan_balance = loan_balance - principal_payment
            if (i == periods .and. abs(loan_balance) <= 100.0_dp * epsilon(amount) * max(1.0_dp, amount)) then
                loan_balance = 0.0_dp
            end if
            result%period(i) = i
            result%interest(i) = round2_r(interest_payment)
            result%principal(i) = round2_r(principal_payment)
            result%balance(i) = round2_r(loan_balance)
        end do
        call set_status(status, cm_success)
    end subroutine loan

end module corpmetrics
