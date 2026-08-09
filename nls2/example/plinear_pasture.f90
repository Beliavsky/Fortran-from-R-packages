module pasture_example_model
    use nls2, only : dp
    implicit none
contains
    subroutine pasture_basis(x, theta, basis, ierr)
        real(dp), intent(in) :: x(:,:), theta(:)
        real(dp), intent(out) :: basis(:,:)
        integer, intent(out) :: ierr
        basis(:,1) = 1.0_dp
        basis(:,2) = -exp(-exp(theta(1) + theta(2) * log(x(:,1))))
        ierr = 0
    end subroutine pasture_basis
end module pasture_example_model

program plinear_pasture
    use nls2
    use pasture_example_model
    implicit none
    integer :: i
    real(dp) :: x(20,1), y(20)
    type(nls_control) :: ctl
    type(nls_result) :: result

    do i = 1, 20
        x(i,1) = 0.5_dp + 0.5_dp * real(i,dp)
        y(i) = 3.0_dp + 70.0_dp * (-exp(-exp(-1.2_dp + 1.5_dp*log(x(i,1)))))
    end do
    ctl%scale_offset = 1.0_dp
    call fit_plinear(pasture_basis, x, y, [-0.5_dp, 1.0_dp], 2, result, ctl)
    print '(a,2f14.8)', 'nonlinear: ', result%par
    print '(a,2f14.8)', 'linear:    ', result%linear_par
    print '(a,es14.6)', 'rss:       ', result%rss
end program plinear_pasture
