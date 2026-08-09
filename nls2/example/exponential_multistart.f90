module exponential_example_model
    use nls2, only : dp
    implicit none
contains
    subroutine exp_model(x, par, yhat, ierr)
        real(dp), intent(in) :: x(:,:), par(:)
        real(dp), intent(out) :: yhat(:)
        integer, intent(out) :: ierr
        yhat = par(1) * exp(par(2) * x(:,1))
        ierr = 0
    end subroutine exp_model
end module exponential_example_model

program exponential_multistart
    use nls2
    use exponential_example_model
    implicit none
    integer :: i
    real(dp) :: x(21,1), y(21), bounds(2,2)
    type(nls_control) :: ctl
    type(nls2_search_result) :: result

    do i = 1, 21
        x(i,1) = real(i-1,dp) / 20.0_dp
        y(i) = 2.0_dp * exp(0.7_dp * x(i,1))
    end do
    bounds(1,:) = [0.5_dp, -1.0_dp]
    bounds(2,:) = [4.0_dp,  2.0_dp]
    ctl%maxiter = 20
    ctl%scale_offset = 1.0_dp
    call seed_rng(1234)
    call nls2_fit(exp_model, x, y, bounds, 'default', result, ctl)
    print '(a,2f14.8)', 'parameters: ', result%best%par
    print '(a,es14.6)', 'rss:        ', result%best%rss
end program exponential_multistart
