program dense_fpca_example
    use fdapace, only : dp, fpca_dense, fpca_result
    implicit none

    real(dp) :: t(21)
    real(dp) :: y(6, 21)
    real(dp) :: score(6)
    type(fpca_result) :: fit
    integer :: i

    do i = 1, size(t)
        t(i) = real(i - 1, dp) / real(size(t) - 1, dp)
    end do
    score = [-1.5_dp, -0.9_dp, -0.3_dp, 0.3_dp, 0.9_dp, 1.5_dp]
    do i = 1, size(y, 1)
        y(i, :) = 2.0_dp + t + score(i) * cos(acos(-1.0_dp) * t)
    end do

    call fpca_dense(y, t, fit, fve_threshold=0.95_dp, max_k=4, assume_error=.false.)
    write (*, '(a,i0)') 'Selected components: ', fit%select_k
    write (*, '(a,f10.6)') 'Explained fraction: ', fit%fve
    write (*, '(a,f10.6)') 'Leading eigenvalue: ', fit%lambda(1)
end program dense_fpca_example
