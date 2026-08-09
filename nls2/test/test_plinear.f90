module test_plinear_models
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
end module test_plinear_models

program test_plinear
    use nls2
    use test_plinear_models
    implicit none
    integer :: i
    real(dp) :: x(20,1), y(20), bounds(2,2)
    type(nls_result) :: r
    type(nls_control) :: ctl
    type(nls2_search_result) :: sr

    do i = 1, 20
        x(i,1) = 0.5_dp + 0.5_dp * real(i,dp)
        y(i) = 3.0_dp + 70.0_dp * (-exp(-exp(-1.2_dp + 1.5_dp*log(x(i,1)))))
    end do
    ctl%scale_offset = 1.0_dp
    call fit_plinear(pasture_basis, x, y, [-0.5_dp, 1.0_dp], 2, r, ctl)
    call assert_true(r%status == nls2_ok .and. r%converged, 'plinear status')
    call assert_close(r%par(1), -1.2_dp, 3.0e-6_dp, 'theta1')
    call assert_close(r%par(2), 1.5_dp, 3.0e-6_dp, 'theta2')
    call assert_close(r%linear_par(1), 3.0_dp, 3.0e-6_dp, 'linear1')
    call assert_close(r%linear_par(2), 70.0_dp, 3.0e-5_dp, 'linear2')

    bounds(1,:) = [-1.2_dp, 1.5_dp]
    bounds(2,:) = [-1.2_dp, 1.5_dp]
    call nls2_fit_plinear(pasture_basis, x, y, bounds, 2, 'plinear-brute-force', sr)
    call assert_true(sr%best%iterations == 0, 'plinear brute does not optimize')
    call assert_true(sr%best%rss < 1.0e-22_dp, 'plinear brute exact')
contains
    subroutine assert_true(ok, label)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: label
        if (.not. ok) error stop 'test_plinear: '//label
    end subroutine assert_true
    subroutine assert_close(x, y, tol, label)
        real(dp), intent(in) :: x, y, tol
        character(len=*), intent(in) :: label
        if (abs(x-y) > tol) error stop 'test_plinear: '//label
    end subroutine assert_close
end program test_plinear
