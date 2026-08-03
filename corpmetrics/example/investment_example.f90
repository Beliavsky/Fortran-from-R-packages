! SPDX-License-Identifier: GPL-2.0-or-later
program investment_example
    use corpmetrics, only : dp, investment_result, idm, cm_success
    implicit none

    type(investment_result) :: result
    integer :: status

    call idm([-1000.0_dp, 300.0_dp, 400.0_dp, 500.0_dp], &
        [0.0_dp, 0.08_dp, 0.09_dp, 0.10_dp], result, status)
    if (status /= cm_success) error stop 'idm failed'

    print '(a,f12.4)', 'NPV: ', result%npv
    print '(a,f12.6)', 'IRR: ', result%irr
end program investment_example
