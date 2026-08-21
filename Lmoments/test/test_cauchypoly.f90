program test_cauchypoly
    use lmoments, only: dp, qcauchypoly, pcauchypoly, dcauchypoly
    implicit none
    real(dp), parameter :: param(4) = [-1.211357448463634_dp, 3.756048230260597_dp, &
        4.533333333333336_dp, 0.258792853488702_dp]
    real(dp), parameter :: probs(5) = [0.1_dp, 0.3_dp, 0.5_dp, 0.8_dp, 0.95_dp]
    real(dp), parameter :: qref(5) = [-1.5869017968250274_dp, 0.1354330066111427_dp, &
        1.7999999999999987_dp, 5.051012273767175_dp, 8.082175474318639_dp]
    real(dp), parameter :: dref(5) = [0.07589102345421236_dp, 0.12956330217635498_dp, &
        0.10986109580510724_dp, 0.07483565851024186_dp, 0.021933547047587418_dp]
    real(dp) :: q, pp, dd
    integer :: i

    do i = 1, size(probs)
        q = qcauchypoly(probs(i), param)
        pp = pcauchypoly(q, param)
        dd = dcauchypoly(q, param)
        call assert_close(q, qref(i), 2.0e-12_dp, 'cauchy polynomial quantile')
        call assert_close(pp, probs(i), 2.0e-14_dp, 'cauchy polynomial inversion')
        call assert_close(dd, dref(i), 3.0e-12_dp, 'cauchy polynomial density')
    end do
    print '(a)', 'test_cauchypoly: PASS'
contains
    subroutine assert_close(a, b, tol, label)
        real(dp), intent(in) :: a, b, tol
        character(*), intent(in) :: label
        if (abs(a-b) > tol) then
            print *, label, a, b
            error stop 1
        end if
    end subroutine
end program test_cauchypoly
