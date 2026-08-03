! SPDX-License-Identifier: GPL-2.0-or-later
program demo_corpmetrics
    use corpmetrics, only : dp, balance_sheet_result, capm_result, ddm_result, &
        fixed_income_result, investment_result, income_statement_result, loan_result, &
        balsh, capm, ddm, fis, idm, insta, loan, valuation_label, ddm_model_label
    implicit none

    type(balance_sheet_result) :: bs
    type(capm_result) :: cp
    type(ddm_result) :: dv
    type(fixed_income_result) :: fi
    type(investment_result) :: inv
    type(income_statement_result) :: inc
    type(loan_result) :: ln
    real(dp) :: rf(5), ri(5), rm(5)
    integer :: status, i

    call balsh(2450000.0_dp, 770000.0_dp, 450000.0_dp, 1180000.0_dp, 490000.0_dp, bs, status)
    print '(a,f12.2)', 'Working capital: ', bs%working_capital
    print '(a,f12.4)', 'Current ratio:   ', bs%current_ratio

    rf = 0.03_dp
    ri = [0.10_dp, 0.14_dp, 0.08_dp, 0.16_dp, 0.12_dp]
    rm = [0.06_dp, 0.11_dp, 0.05_dp, 0.13_dp, 0.09_dp]
    call capm(rf, ri, rm, cp, status)
    print '(a,f10.5)', 'CAPM required return: ', cp%required_return
    print '(a,a)', 'Valuation: ', valuation_label(cp%valuation)

    call ddm(2.0_dp, 0.12_dp, dv, g1=0.08_dp, g2=0.04_dp, period=3, status=status)
    print '(a,a)', 'DDM model: ', ddm_model_label(dv%model)
    print '(a,f12.4)', 'Stock value: ', dv%value

    call fis(1000.0_dp, 0.08_dp, 0.12_dp, 2, fi, semiannual=.true., status=status)
    print '(a,f12.4)', 'Bond price: ', fi%price
    print '(a,f12.6)', 'Macaulay duration: ', fi%macaulay_duration

    call idm([-100.0_dp, 120.0_dp], [0.0_dp, 0.1_dp], inv, status)
    print '(a,f12.4)', 'NPV: ', inv%npv
    print '(a,f12.4,a)', 'IRR: ', inv%irr_percent, '%'

    call insta(25000000.0_dp, 19850000.0_dp, 1000000.0_dp, inc, &
        preferred_dividend=100000.0_dp, shares=100000.0_dp, price_per_share=120.0_dp, status=status)
    print '(a,f10.4)', 'Gross margin (%): ', inc%gross_profit_margin_percent
    print '(a,f10.4)', 'EPS: ', inc%earnings_per_share
    print '(a,f10.4)', 'P/E: ', inc%price_earnings_ratio

    call loan(1000.0_dp, 0.20_dp, 3, ln, status)
    print '(a,f12.4)', 'Loan installment: ', ln%installment
    print '(a)', 'Period     Interest    Principal      Balance'
    do i = 1, size(ln%period)
        print '(i6,3f13.2)', ln%period(i), ln%interest(i), ln%principal(i), ln%balance(i)
    end do
end program demo_corpmetrics
