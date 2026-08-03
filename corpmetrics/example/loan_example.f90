! SPDX-License-Identifier: GPL-2.0-or-later
program loan_example
    use corpmetrics, only : dp, loan_result, loan, cm_success
    implicit none

    type(loan_result) :: result
    integer :: status, i

    call loan(100000.0_dp, 0.05_dp, 4, result, status)
    if (status /= cm_success) error stop 'loan failed'

    print '(a,f14.2)', 'Installment: ', result%installment
    print '(a,f14.2)', 'Total repayment: ', result%total_repayment
    do i = 1, size(result%period)
        print '(i4,3f14.2)', result%period(i), result%interest(i), &
            result%principal(i), result%balance(i)
    end do
end program loan_example
