program test_monotonicity
    use goftest, only : dp, p_cvm, p_ad
    implicit none

    integer :: i, j
    integer, parameter :: nvals(7) = [1, 2, 5, 10, 50, 100, 500]
    real(dp) :: x, prev, cur

    prev = 0.0_dp
    do i = 1, 120
        x = 0.005_dp + real(i - 1, dp) * 0.025_dp
        cur = p_cvm(x)
        if (cur + 2.0e-10_dp < prev) error stop 1
        prev = cur
    end do

    prev = 0.0_dp
    do i = 1, 120
        x = 0.02_dp + real(i - 1, dp) * 0.05_dp
        cur = p_ad(x)
        if (cur + 1.0e-12_dp < prev) error stop 1
        prev = cur
    end do

    do j = 1, size(nvals)
        prev = 0.0_dp
        do i = 1, 160
            x = 0.005_dp + real(i - 1, dp) * 0.05_dp
            cur = p_cvm(x, n=nvals(j))
            if (cur + 2.0e-9_dp < prev) error stop 1
            prev = cur
        end do
    end do

    print '(a)', 'test_monotonicity: PASS'
end program test_monotonicity
