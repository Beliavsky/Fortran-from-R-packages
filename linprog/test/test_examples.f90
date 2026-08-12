program test_examples
    use linprog
    implicit none
    real(dp) :: c1(3), b1(3), a1(3,3), c2(2), b2(3), a2(3,2)
    type(linprog_result) :: r
    type(linprog_control) :: ctl

    c1 = [1800.0_dp, 600.0_dp, 600.0_dp]
    b1 = [40.0_dp, 90.0_dp, 2500.0_dp]
    a1 = reshape([0.7_dp,1.5_dp,50.0_dp, 0.35_dp,1.0_dp,12.5_dp, &
                  0.0_dp,3.0_dp,20.0_dp], [3,3])
    ctl%maximum = .true.
    call solveLP(c1, b1, a1, r, ctl)
    call assert_true(r%status == 0, 'production status')
    call assert_close(r%opt, 93600.0_dp, 1.0e-8_dp, 'production objective')
    call assert_vec(r%solution, [44.0_dp,24.0_dp,0.0_dp], 1.0e-8_dp, 'production solution')
    call assert_vec(r%con_dual, [0.0_dp,240.0_dp,28.8_dp], 1.0e-8_dp, 'production duals')
    call assert_true(r%iter1 == 0 .and. r%iter2 == 2, 'production iterations')
    call assert_close(r%allvar_max_c(3), 1296.0_dp, 1.0e-8_dp, 'Pigs max c')

    c2 = [2.5_dp, 2.0_dp]
    b2 = [-10.0_dp, -1.5_dp, 12.0_dp]
    a2 = reshape([-1.6_dp,-0.5_dp,2.0_dp, -2.4_dp,-0.2_dp,2.0_dp], [3,2])
    ctl = linprog_control()
    call solveLP(c2, b2, a2, r, ctl)
    call assert_true(r%status == 0, 'feed status')
    call assert_close(r%opt, 10.454545455_dp, 1.0e-9_dp, 'feed objective')
    call assert_vec(r%solution, [1.818181818_dp,2.954545455_dp], 1.0e-9_dp, 'feed solution')
    call assert_true(r%iter1 == 2 .and. r%iter2 == 0, 'feed iterations')
    call assert_vec(r%con_dual, [0.568181818181818_dp,3.18181818181818_dp,0.0_dp], &
        1.0e-10_dp, 'feed duals')
    print *, 'test_examples: PASS'
contains
    subroutine assert_true(ok, msg)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: msg
        if (.not. ok) error stop msg
    end subroutine
    subroutine assert_close(x, y, tol, msg)
        real(dp), intent(in) :: x, y, tol
        character(len=*), intent(in) :: msg
        if (abs(x-y) > tol) then
            print *, trim(msg), x, y
            error stop 'assert_close failed'
        end if
    end subroutine
    subroutine assert_vec(x, y, tol, msg)
        real(dp), intent(in) :: x(:), y(:), tol
        character(len=*), intent(in) :: msg
        if (size(x) /= size(y) .or. maxval(abs(x-y)) > tol) then
            print *, trim(msg), x, y
            error stop 'assert_vec failed'
        end if
    end subroutine
end program test_examples
