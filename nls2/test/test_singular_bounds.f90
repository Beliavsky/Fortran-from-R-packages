module test_sb_models
    use nls2, only : dp
    implicit none
contains
    subroutine singular_model(x, par, yhat, ierr)
        real(dp), intent(in) :: x(:,:), par(:)
        real(dp), intent(out) :: yhat(:)
        integer, intent(out) :: ierr
        yhat = (par(1) + 2.0_dp*par(2)) * x(:,1)
        ierr = 0
    end subroutine singular_model
    subroutine onepar_model(x, par, yhat, ierr)
        real(dp), intent(in) :: x(:,:), par(:)
        real(dp), intent(out) :: yhat(:)
        integer, intent(out) :: ierr
        yhat = par(1) * x(:,1)
        ierr = 0
    end subroutine onepar_model
end module test_sb_models

program test_singular_bounds
    use nls2
    use test_sb_models
    implicit none
    integer :: i
    real(dp) :: x(9,1), y(9), one_start(1,2), bounds(2,1), y2(9)
    type(nls_result) :: r
    type(nls2_search_result) :: sr

    do i = 1, 9
        x(i,1) = 1.0_dp
        y(i) = real(i,dp)
        y2(i) = 5.0_dp
    end do
    call fit_nls(singular_model, x, y, [1.0_dp, 1.0_dp], r)
    call assert_true(r%status == nls2_singular, 'singular gradient detected')
    one_start(1,:) = [1.0_dp, 1.0_dp]
    call nls2_fit(singular_model, x, y, one_start, 'brute-force', sr)
    call assert_true(sr%status == nls2_ok, 'brute allows singular Jacobian point')
    call assert_true(sr%best%iterations == 0, 'singular brute evaluation only')

    bounds(:,1) = [0.0_dp, 2.0_dp]
    call nls2_fit(onepar_model, x, y2, bounds, 'port', sr, lower=[0.0_dp], upper=[2.0_dp])
    call assert_true(sr%best%par(1) <= 2.0_dp + 1.0e-12_dp, 'port compatibility upper bound')
    call assert_close(sr%best%par(1), 2.0_dp, 1.0e-9_dp, 'bounded optimum')
contains
    subroutine assert_true(ok, label)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: label
        if (.not. ok) error stop 'test_singular_bounds: '//label
    end subroutine assert_true
    subroutine assert_close(x, y, tol, label)
        real(dp), intent(in) :: x, y, tol
        character(len=*), intent(in) :: label
        if (abs(x-y) > tol) error stop 'test_singular_bounds: '//label
    end subroutine assert_close
end program test_singular_bounds
