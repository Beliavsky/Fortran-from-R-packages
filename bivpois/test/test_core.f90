program test_core
    use bivpois, only : dp, dbp_scalar, dbp, bp_grid, bp_table, &
                        bp_probability_grid, make_bp_table
    implicit none
    real(dp), parameter :: lambda(3) = [3.0_dp, 5.0_dp, 2.0_dp]
    integer, parameter :: x(6) = [0, 1, 0, 1, 3, 7]
    integer, parameter :: y(6) = [0, 0, 1, 1, 4, 8]
    real(dp), parameter :: ref(6) = [ &
        4.5399929762484854e-5_dp, 1.3619978928745447e-4_dp, &
        2.2699964881242410e-4_dp, 7.717988059622419e-4_dp, &
        1.7540452030110038e-2_dp, 1.5765439027418143e-2_dp ]
    real(dp) :: p(6), total, indep, expected
    real(dp) :: lambda0(3)
    integer :: i, j
    type(bp_grid) :: grid
    type(bp_table) :: tab

    call dbp(x, y, lambda, p, logged=.false.)
    call assert_close_vec(p, ref, 2.0e-14_dp, "reference probabilities")

    total = 0.0_dp
    do j = 0, 40
        do i = 0, 40
            total = total + dbp_scalar(i, j, lambda, .false.)
        end do
    end do
    call assert_close(total, 1.0_dp, 2.0e-12_dp, "probability mass")

    lambda0 = [2.0_dp, 4.0_dp, 0.0_dp]
    indep = dbp_scalar(3, 5, lambda0, .false.)
    expected = exp(-2.0_dp) * 2.0_dp**3 / 6.0_dp * &
               exp(-4.0_dp) * 4.0_dp**5 / 120.0_dp
    call assert_close(indep, expected, 1.0e-15_dp, "independence reduction")

    grid = bp_probability_grid([1, 2, 3], [2, 3, 4], lambda, padding=1)
    if (lbound(grid%probability, 1) /= 0 .or. ubound(grid%probability, 1) /= 4) &
        error stop "grid x bounds"
    if (lbound(grid%probability, 2) /= 1 .or. ubound(grid%probability, 2) /= 5) &
        error stop "grid y bounds"

    tab = make_bp_table([0, 1, 1, 2], [1, 1, 2, 2])
    if (sum(tab%count) /= 4) error stop "table total"
    if (tab%count(1, 1) /= 1 .or. tab%count(1, 2) /= 1) error stop "table cells"

    print '(a)', 'test_core: PASS'

contains

    subroutine assert_close(a, b, tol, label)
        real(dp), intent(in) :: a, b, tol
        character(*), intent(in) :: label
        if (abs(a - b) > tol) then
            print '(a,2es24.15)', trim(label)//": ", a, b
            error stop 1
        end if
    end subroutine assert_close

    subroutine assert_close_vec(a, b, tol, label)
        real(dp), intent(in) :: a(:), b(:), tol
        character(*), intent(in) :: label
        if (size(a) /= size(b)) error stop "assert_close_vec: size mismatch"
        if (maxval(abs(a - b)) > tol) then
            print '(a,es24.15)', trim(label)//" max error: ", maxval(abs(a - b))
            error stop 1
        end if
    end subroutine assert_close_vec

end program test_core
