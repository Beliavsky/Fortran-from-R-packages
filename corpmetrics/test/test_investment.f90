! SPDX-License-Identifier: GPL-2.0-or-later
program test_investment
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use corpmetrics, only : dp, cm_success, cm_dimension_mismatch, cm_root_not_bracketed, &
        investment_result, idm, net_present_value, internal_rate_of_return
    implicit none

    type(investment_result) :: result
    real(dp) :: irr
    integer :: status

    call idm([-100.0_dp, 120.0_dp], [0.0_dp, 0.1_dp], result, status)
    call assert_int(status, cm_success, 'simple idm status')
    call assert_close(result%npv, 9.09090909090909_dp, 1.0e-12_dp, 'simple npv')
    call assert_close(result%irr, 0.2_dp, 1.0e-9_dp, 'simple irr')
    call assert_close(result%irr_percent, 20.0_dp, 1.0e-9_dp, 'simple irr percent')

    call idm([-1000.0_dp, 300.0_dp, 400.0_dp, 500.0_dp], &
        [0.0_dp, 0.08_dp, 0.09_dp, 0.10_dp], result, status)
    call assert_int(status, cm_success, 'multi-period idm status')
    call assert_close(result%npv, -9.89282446480945_dp, 1.0e-12_dp, 'multi-period npv')
    call assert_close(result%irr, 0.08896339469329151_dp, 2.0e-9_dp, 'multi-period irr')

    call assert_close(net_present_value([-100.0_dp, 120.0_dp], [0.0_dp, 0.1_dp]), &
        9.09090909090909_dp, 1.0e-12_dp, 'public npv')

    call internal_rate_of_return([100.0_dp, 20.0_dp], irr, status)
    call assert_int(status, cm_root_not_bracketed, 'unbracketed irr status')
    if (.not. ieee_is_nan(irr)) error stop 'unbracketed irr should be NaN'

    call idm([-100.0_dp, 120.0_dp], [0.0_dp], result, status)
    call assert_int(status, cm_dimension_mismatch, 'idm dimension status')

    print '(a)', 'test_investment: PASS'

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

end program test_investment
