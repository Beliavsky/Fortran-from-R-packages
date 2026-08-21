program test_normpoly
    use lmoments, only: dp, lmom2normpoly4, lmom2normpoly6, qnormpoly, pnormpoly, &
        dnormpoly, data2normpoly4
    implicit none
    real(dp), parameter :: lm(6) = [2.27_dp, 2.081111111111111_dp, 0.243333333333333_dp, &
        0.398571428571429_dp, 0.092222222222222_dp, 0.25_dp]
    real(dp), parameter :: p4ref(4) = [6.996167879087116_dp, -14.319002424840892_dp, &
        7.299999999999990_dp, 5.762149715989586_dp]
    real(dp), parameter :: p6ref(6) = [17.375122558869048_dp, -71.25739252215854_dp, &
        145.72144221326127_dp, -158.68096147550736_dp, 58.09999999999986_dp, &
        10.148917818159129_dp]
    real(dp), parameter :: probs(5) = [0.1_dp, 0.3_dp, 0.5_dp, 0.8_dp, 0.95_dp]
    real(dp), parameter :: qref(5) = [-1.7472243528258025_dp, 0.3357928862694166_dp, &
        1.661666666666667_dp, 5.062513491217368_dp, 9.459258434871124_dp]
    real(dp), parameter :: dref(5) = [0.050064905252669815_dp, 0.15074912989475633_dp, &
        0.13468802593066165_dp, 0.05573234222084805_dp, 0.018043824208888632_dp]
    real(dp), parameter :: xdata(10) = [3.2_dp, -1.1_dp, 0.4_dp, 5.7_dp, 2.2_dp, &
        1.8_dp, 9.1_dp, -3.0_dp, 4.4_dp, 0.0_dp]
    real(dp) :: p4(4), p6(6), q, pp, dd
    integer :: i

    p4 = lmom2normpoly4(lm)
    p6 = lmom2normpoly6(lm)
    call assert_vec(p4, p4ref, 2.0e-12_dp, 'normpoly4 parameters')
    call assert_vec(p6, p6ref, 2.0e-11_dp, 'normpoly6 parameters')
    call assert_vec(data2normpoly4(xdata), p4ref, 2.0e-12_dp, 'data2normpoly4')

    do i = 1, size(probs)
        q = qnormpoly(probs(i), p4)
        pp = pnormpoly(q, p4)
        dd = dnormpoly(q, p4)
        call assert_close(q, qref(i), 2.0e-11_dp, 'normal polynomial quantile')
        call assert_close(pp, probs(i), 2.0e-14_dp, 'normal polynomial inversion')
        call assert_close(dd, dref(i), 3.0e-11_dp, 'normal polynomial density')
    end do
    print '(a)', 'test_normpoly: PASS'
contains
    subroutine assert_vec(a, b, tol, label)
        real(dp), intent(in) :: a(:), b(:), tol
        character(*), intent(in) :: label
        if (maxval(abs(a-b)) > tol) then
            print *, label, a, b
            error stop 1
        end if
    end subroutine
    subroutine assert_close(a, b, tol, label)
        real(dp), intent(in) :: a, b, tol
        character(*), intent(in) :: label
        if (abs(a-b) > tol) then
            print *, label, a, b
            error stop 1
        end if
    end subroutine
end program test_normpoly
