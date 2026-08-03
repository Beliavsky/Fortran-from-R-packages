! SPDX-License-Identifier: GPL-2.0-or-later
program test_loan
    use corpmetrics, only : dp, cm_success, loan_result, loan, round2_r
    implicit none

    type(loan_result) :: result
    integer :: status

    call loan(1000.0_dp, 0.20_dp, 3, result, status)
    call assert_int(status, cm_success, 'loan status')
    call assert_close(result%installment, 474.72527472527474_dp, 1.0e-12_dp, 'installment')
    call assert_close(result%total_repayment, 1424.1758241758243_dp, 1.0e-12_dp, 'total')
    call assert_array(result%interest, [200.0_dp, 145.05_dp, 79.12_dp], 'interest')
    call assert_array(result%principal, [274.73_dp, 329.67_dp, 395.60_dp], 'principal')
    call assert_array(result%balance, [725.27_dp, 395.60_dp, 0.0_dp], 'balance')

    call loan(1200.0_dp, 0.0_dp, 12, result, status)
    call assert_int(status, cm_success, 'zero-rate status')
    call assert_close(result%installment, 100.0_dp, 1.0e-13_dp, 'zero-rate installment')
    call assert_close(result%balance(12), 0.0_dp, 1.0e-13_dp, 'zero-rate final balance')

    call assert_close(round2_r(2.125_dp), 2.12_dp, 1.0e-13_dp, 'round half even')
    call assert_close(round2_r(2.375_dp), 2.38_dp, 1.0e-13_dp, 'round half odd')

    print '(a)', 'test_loan: PASS'

contains

    subroutine assert_array(actual, expected, label)
        real(dp), intent(in) :: actual(:), expected(:)
        character(len=*), intent(in) :: label

        if (size(actual) /= size(expected)) error stop trim(label) // ': size mismatch'
        if (any(abs(actual - expected) > 1.0e-12_dp)) then
            print '(a)', trim(label) // ': array mismatch'
            error stop 1
        end if
    end subroutine assert_array

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

end program test_loan
