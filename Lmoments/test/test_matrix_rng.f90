program test_matrix_rng
    use lmoments, only: dp, lmoments_matrix, rnormpoly, rcauchypoly, pnormpoly, pcauchypoly
    implicit none
    real(dp) :: data(10,2), lm(2,4), rn(2000), rc(2000)
    real(dp), parameter :: np(4) = [0.0_dp, 1.0_dp, 0.2_dp, 0.8_dp]
    real(dp), parameter :: cp(4) = [0.0_dp, 1.0_dp, 0.2_dp, 0.2_dp]
    integer :: i, info, inside_n, inside_c

    do i = 1, 10
        data(i,1) = real(i,dp)
        data(i,2) = 2.0_dp * real(i,dp) - 3.0_dp
    end do
    call lmoments_matrix(data, lm, info)
    call assert_true(info == 0, 'matrix info')
    call assert_close(lm(2,1), 2.0_dp * lm(1,1) - 3.0_dp, 1.0e-13_dp, 'location equivariance')
    call assert_close(lm(2,2), 2.0_dp * lm(1,2), 1.0e-13_dp, 'scale equivariance')

    call random_seed()
    call rnormpoly(np, rn)
    call rcauchypoly(cp, rc)
    inside_n = 0
    inside_c = 0
    do i = 1, size(rn)
        if (pnormpoly(rn(i), np) > 0.01_dp .and. pnormpoly(rn(i), np) < 0.99_dp) inside_n = inside_n + 1
        if (pcauchypoly(rc(i), cp) > 0.01_dp .and. pcauchypoly(rc(i), cp) < 0.99_dp) inside_c = inside_c + 1
    end do
    call assert_true(inside_n > 1900, 'normal PIT range')
    call assert_true(inside_c > 1900, 'cauchy PIT range')

    print '(a)', 'test_matrix_rng: PASS'
contains
    subroutine assert_close(a, b, tol, label)
        real(dp), intent(in) :: a, b, tol
        character(*), intent(in) :: label
        if (abs(a-b) > tol) then
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
end program test_matrix_rng
