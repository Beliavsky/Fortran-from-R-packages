module test_search_models
    use nls2, only : dp
    implicit none
contains
    subroutine linear_model(x, par, yhat, ierr)
        real(dp), intent(in) :: x(:,:), par(:)
        real(dp), intent(out) :: yhat(:)
        integer, intent(out) :: ierr
        yhat = par(1) + par(2) * x(:,1)
        ierr = 0
    end subroutine linear_model
end module test_search_models

program test_search
    use nls2
    use test_search_models
    implicit none
    integer :: i
    real(dp) :: x(9,1), y(9), bounds(2,2), rows(3,2)
    real(dp), allocatable :: grid(:,:)
    type(nls2_search_result) :: r
    type(nls_control) :: ctl

    do i = 1, 9
        x(i,1) = real(i-5,dp)
        y(i) = 1.0_dp - x(i,1) / 3.0_dp
    end do
    bounds(1,:) = [-1.0_dp, -1.0_dp]
    bounds(2,:) = [ 1.0_dp,  1.0_dp]
    call make_grid(bounds(1,:), bounds(2,:), 10, grid)
    call assert_true(size(grid,1) == 16, 'source-compatible ceil grid size')
    call assert_close(grid(1,1), -1.0_dp, 0.0_dp, 'grid first')
    call assert_close(grid(2,1), -1.0_dp/3.0_dp, 1.0e-15_dp, 'first coordinate varies fastest')

    ctl%maxiter = 10
    call nls2_fit(linear_model, x, y, bounds, 'brute-force', r, ctl)
    call assert_true(r%n_candidates == 16, 'brute candidate count')
    call assert_true(r%best%iterations == 0, 'brute does not optimize')
    call assert_true(r%best%rss < 1.0e-26_dp, 'brute exact candidate')

    ctl%maxiter = 12
    ctl%tol = 1.0e-7_dp
    call seed_rng(12345)
    call nls2_fit(linear_model, x, y, bounds, 'default', r, ctl)
    call assert_true(r%n_candidates == 12, 'default bounds random starts')
    call assert_true(r%best%rss < 1.0e-15_dp, 'default multi-start optimizes')

    rows(1,:) = [-0.5_dp, 0.5_dp]
    rows(2,:) = [ 0.0_dp, 0.0_dp]
    rows(3,:) = [ 1.0_dp,-1.0_dp]
    ctl%maxiter = 2
    call nls2_fit(linear_model, x, y, rows, 'random-search', r, ctl)
    call assert_true(r%n_candidates == 3, 'current source uses all explicit rows')
contains
    subroutine assert_true(ok, label)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: label
        if (.not. ok) error stop 'test_search: '//label
    end subroutine assert_true
    subroutine assert_close(x, y, tol, label)
        real(dp), intent(in) :: x, y, tol
        character(len=*), intent(in) :: label
        if (abs(x-y) > tol) error stop 'test_search: '//label
    end subroutine assert_close
end program test_search
