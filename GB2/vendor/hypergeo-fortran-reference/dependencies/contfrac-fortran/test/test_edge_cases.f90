! SPDX-License-Identifier: GPL-2.0-only
program test_edge_cases
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_positive_inf
    use contfrac, only : dp, cf, gcf, as_cf, contfrac_info
    implicit none

    real(dp) :: value, inf
    complex(dp) :: cvalue, cexpected
    real(dp), allocatable :: a(:), terms(:)
    type(contfrac_info) :: info
    integer :: n_used

    inf = ieee_value(0.0_dp, ieee_positive_inf)
    a = [1.0_dp, 2.0_dp, inf, 99.0_dp]
    value = cf(a)
    call assert_close(value, 1.5_dp, 2.0e-15_dp, "infinite coefficient truncation")

    value = gcf([2.0_dp], [3.0_dp], b0=1.0_dp, finite=.false., info=info)
    if (.not. ieee_is_nan(value)) error stop "nonconverged GCF should return NaN"
    if (info%converged) error stop "nonconverged GCF status incorrect"

    value = gcf([2.0_dp], [3.0_dp], b0=1.0_dp, finite=.false., tol=1.0e30_dp, info=info)
    call assert_close(value, 5.0_dp / 3.0_dp, 5.0e-15_dp, "loose acceptance tolerance")


    cexpected = cmplx(1.0_dp, 0.5_dp, kind=dp) + &
        1.0_dp / cmplx(2.0_dp, -0.25_dp, kind=dp)
    cvalue = cf([cmplx(1.0_dp, 0.5_dp, kind=dp), &
                 cmplx(2.0_dp, -0.25_dp, kind=dp)], finite=.true.)
    if (abs(cvalue - cexpected) > 5.0e-15_dp) error stop "complex CF wrapper incorrect"

    terms = as_cf(1.5_dp, 6, n_used)
    if (n_used /= 2) error stop "rational as_cf n_used incorrect"
    call assert_close(terms(1), 1.0_dp, 0.0_dp, "rational term 1")
    call assert_close(terms(2), 2.0_dp, 0.0_dp, "rational term 2")
    if (.not. ieee_is_nan(terms(3))) error stop "rational tail should be NaN"

    print '(a)', "test_edge_cases: PASS"

contains

    subroutine assert_close(actual, expected, atol, label)
        real(dp), intent(in) :: actual, expected, atol
        character(len=*), intent(in) :: label
        if (abs(actual - expected) > atol) then
            print '(a,2es25.16)', trim(label) // ": ", actual, expected
            error stop 1
        end if
    end subroutine assert_close

end program test_edge_cases
