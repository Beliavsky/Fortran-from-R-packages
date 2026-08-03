! SPDX-License-Identifier: GPL-2.0-or-later
program test_metrics
    use corpmetrics, only : dp, cm_success, capm_undervalued, ddm_zero_growth, &
        ddm_gordon, ddm_differential, balance_sheet_result, capm_result, ddm_result, &
        fixed_income_result, income_statement_result, balsh, capm, ddm, fis, insta
    implicit none

    type(balance_sheet_result) :: bs
    type(capm_result) :: cp
    type(ddm_result) :: dv
    type(fixed_income_result) :: fi
    type(income_statement_result) :: inc
    real(dp) :: rf(5), ri(5), rm(5)
    integer :: status

    call balsh(2450000.0_dp, 770000.0_dp, 450000.0_dp, 1180000.0_dp, 490000.0_dp, bs, status)
    call assert_int(status, cm_success, 'balsh status')
    call assert_close(bs%working_capital, 280000.0_dp, 1.0e-12_dp, 'working capital')
    call assert_close(bs%current_ratio, 1.5714285714285714_dp, 1.0e-13_dp, 'current ratio')
    call assert_close(bs%acid_test_ratio, 0.6530612244897959_dp, 1.0e-13_dp, 'acid ratio')
    call assert_close(bs%leverage_ratio, 0.5186335403726708_dp, 1.0e-13_dp, 'leverage')
    call assert_close(bs%debt_to_equity_ratio, 1.0774193548387097_dp, 1.0e-13_dp, 'debt/equity')

    rf = 0.03_dp
    ri = [0.10_dp, 0.14_dp, 0.08_dp, 0.16_dp, 0.12_dp]
    rm = [0.06_dp, 0.11_dp, 0.05_dp, 0.13_dp, 0.09_dp]
    call capm(rf, ri, rm, cp, status)
    call assert_int(status, cm_success, 'capm status')
    call assert_close(cp%beta, 0.9375_dp, 1.0e-13_dp, 'beta')
    call assert_close(cp%required_return, 0.084375_dp, 1.0e-13_dp, 'required return')
    call assert_int(cp%valuation, capm_undervalued, 'valuation')

    call ddm(3.0_dp, 0.08_dp, dv, status=status)
    call assert_int(status, cm_success, 'zero-growth status')
    call assert_int(dv%model, ddm_zero_growth, 'zero-growth model')
    call assert_close(dv%value, 37.5_dp, 1.0e-13_dp, 'zero-growth value')

    call ddm(0.8_dp, 0.10_dp, dv, g1=0.04_dp, status=status)
    call assert_int(status, cm_success, 'gordon status')
    call assert_int(dv%model, ddm_gordon, 'gordon model')
    call assert_close(dv%value, 13.866666666666667_dp, 1.0e-13_dp, 'gordon value')

    call ddm(2.0_dp, 0.12_dp, dv, g1=0.08_dp, g2=0.04_dp, period=3, status=status)
    call assert_int(status, cm_success, 'differential status')
    call assert_int(dv%model, ddm_differential, 'differential model')
    call assert_close(dv%value, 29.79077077259476_dp, 1.0e-12_dp, 'differential value')

    call fis(1000.0_dp, 0.08_dp, 0.08_dp, 6, fi, status=status)
    call assert_int(status, cm_success, 'annual fis status')
    call assert_close(fi%price, 1000.0_dp, 1.0e-11_dp, 'annual bond price')
    call assert_close(fi%macaulay_duration, 4.992710037078084_dp, 1.0e-12_dp, 'annual duration')
    call assert_close(fi%modified_duration, 4.622879663961189_dp, 1.0e-12_dp, 'annual modified')

    call fis(1000.0_dp, 0.08_dp, 0.12_dp, 2, fi, semiannual=.true., status=status)
    call assert_int(status, cm_success, 'semiannual fis status')
    call assert_close(fi%price, 930.6978877460067_dp, 1.0e-11_dp, 'semiannual price')
    call assert_close(fi%macaulay_duration, 1.8828878648149816_dp, 1.0e-12_dp, 'semiannual duration')
    call assert_close(fi%modified_duration, 1.6811498792990907_dp, 1.0e-12_dp, 'semiannual modified')

    call insta(25000000.0_dp, 19850000.0_dp, 1000000.0_dp, inc, &
        preferred_dividend=100000.0_dp, shares=100000.0_dp, price_per_share=120.0_dp, status=status)
    call assert_int(status, cm_success, 'insta status')
    call assert_close(inc%gross_profit_margin_percent, 20.6_dp, 1.0e-13_dp, 'gross margin')
    call assert_close(inc%net_profit_margin_percent, 4.0_dp, 1.0e-13_dp, 'net margin')
    call assert_close(inc%earnings_per_share, 9.0_dp, 1.0e-13_dp, 'eps')
    call assert_close(inc%price_earnings_ratio, 13.333333333333334_dp, 1.0e-13_dp, 'pe')

    print '(a)', 'test_metrics: PASS'

contains

    subroutine assert_close(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: label

        if (abs(actual - expected) > tolerance * max(1.0_dp, abs(expected))) then
            print '(a,2es24.15)', trim(label) // ': ', actual, expected
            error stop 1
        end if
    end subroutine assert_close

    subroutine assert_int(actual, expected, label)
        integer, intent(in) :: actual, expected
        character(len=*), intent(in) :: label

        if (actual /= expected) then
            print '(a,2i12)', trim(label) // ': ', actual, expected
            error stop 1
        end if
    end subroutine assert_int

end program test_metrics
