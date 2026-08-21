! SPDX-License-Identifier: GPL-2.0-or-later
program test_stats
    use adequacy_model, only: dp, descriptive_result, descriptive, ttt_curve
    implicit none
    real(dp) :: x(6), y(3)
    real(dp), allocatable :: r(:), t(:)
    type(descriptive_result) :: d

    x = [1.0_dp, 2.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 6.0_dp]
    d = descriptive(x)
    call assert_close(d%mean, 3.0_dp, 1.0e-12_dp, 'mean')
    call assert_close(d%median, 2.5_dp, 1.0e-12_dp, 'median')
    call assert_close(d%variance, 3.2_dp, 1.0e-12_dp, 'variance')
    if (size(d%mode) /= 1) error stop 'mode size'
    call assert_close(d%mode(1), 2.0_dp, 1.0e-12_dp, 'mode')
    if (d%n /= 6) error stop 'n'

    y = [1.0_dp, 2.0_dp, 3.0_dp]
    call ttt_curve(y, r, t)
    call assert_close(r(1), 1.0_dp/3.0_dp, 1.0e-12_dp, 'r1')
    call assert_close(r(3), 1.0_dp, 1.0e-12_dp, 'r3')
    call assert_close(t(1), 0.5_dp, 1.0e-12_dp, 't1')
    call assert_close(t(2), 5.0_dp/6.0_dp, 1.0e-12_dp, 't2')
    call assert_close(t(3), 1.0_dp, 1.0e-12_dp, 't3')

    print '(a)', 'test_stats: PASS'
contains
    subroutine assert_close(a, b, tol, label)
        real(dp), intent(in) :: a, b, tol
        character(len=*), intent(in) :: label
        if (abs(a-b) > tol) then
            print *, trim(label), a, b
            error stop 1
        end if
    end subroutine assert_close
end program test_stats
