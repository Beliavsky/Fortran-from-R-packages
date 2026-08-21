program test_t1
    use lmoments, only: dp, t1_lmoments, t1lmom2cauchypoly4, data2cauchypoly4
    implicit none
    real(dp), parameter :: x(10) = [3.2_dp, -1.1_dp, 0.4_dp, 5.7_dp, 2.2_dp, &
        1.8_dp, 9.1_dp, -3.0_dp, 4.4_dp, 0.0_dp]
    real(dp), parameter :: ref(4) = [2.026666666666666_dp, 1.009523809523810_dp, &
        0.071957671957672_dp, 0.061904761904762_dp]
    real(dp), parameter :: pref(4) = [-1.211357448463634_dp, 3.756048230260597_dp, &
        4.533333333333336_dp, 0.258792853488702_dp]
    real(dp) :: t(4), p(4), p2(4)
    integer :: info

    call t1_lmoments(x, t, info)
    call assert_true(info == 0, 't1 info')
    call assert_vec(t, ref, 3.0e-13_dp, 't1 reference')
    p = t1lmom2cauchypoly4(t)
    p2 = data2cauchypoly4(x)
    call assert_vec(p, pref, 4.0e-12_dp, 'cauchy parameters')
    call assert_vec(p2, pref, 4.0e-12_dp, 'data2 cauchy parameters')
    print '(a)', 'test_t1: PASS'
contains
    subroutine assert_vec(a, b, tol, label)
        real(dp), intent(in) :: a(:), b(:), tol
        character(*), intent(in) :: label
        if (maxval(abs(a - b)) > tol) then
            print *, label, a, b
            error stop 1
        end if
    end subroutine
    subroutine assert_true(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not. ok) then
            print *, label
            error stop 1
        end if
    end subroutine
end program test_t1
