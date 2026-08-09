module test_nls_models
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

    subroutine exp_jac(x, par, jac, ierr)
        real(dp), intent(in) :: x(:,:), par(:)
        real(dp), intent(out) :: jac(:,:)
        integer, intent(out) :: ierr
        jac(:,1) = exp(par(2) * x(:,1))
        jac(:,2) = par(1) * x(:,1) * exp(par(2) * x(:,1))
        ierr = 0
    end subroutine exp_jac
end module test_nls_models

program test_nls
    use nls2
    use test_nls_models
    implicit none
    integer :: i
    real(dp) :: x(21,1), y(21), w(21), ll
    type(nls_result) :: a, b

    do i = 1, 21
        x(i,1) = real(i-1,dp) / 20.0_dp
        y(i) = 2.0_dp * exp(0.7_dp * x(i,1))
        w(i) = 0.5_dp + real(i,dp) / 21.0_dp
    end do
    call fit_nls(exp_model, x, y, [1.0_dp, 0.1_dp], a)
    call assert_true(a%status == nls2_ok .and. a%converged, 'numerical fit status')
    call assert_close(a%par(1), 2.0_dp, 2.0e-6_dp, 'a')
    call assert_close(a%par(2), 0.7_dp, 2.0e-6_dp, 'b')
    call assert_true(a%rss < 1.0e-15_dp, 'numerical rss')

    call fit_nls(exp_model, x, y, [1.2_dp, 0.2_dp], b, weights=w, jacobian=exp_jac)
    call assert_true(b%status == nls2_ok .and. b%converged, 'analytic fit status')
    call assert_close(b%par(1), 2.0_dp, 2.0e-7_dp, 'weighted a')
    call assert_close(b%par(2), 0.7_dp, 2.0e-7_dp, 'weighted b')
    ll = nls_loglik(b, w)
    call assert_true(ll > 0.0_dp, 'exact-fit logLik positive/infinite')
    call assert_true(nls_df_residual(b, w) == 19, 'residual df')
contains
    subroutine assert_true(ok, label)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: label
        if (.not. ok) error stop 'test_nls: '//label
    end subroutine assert_true
    subroutine assert_close(x, y, tol, label)
        real(dp), intent(in) :: x, y, tol
        character(len=*), intent(in) :: label
        if (abs(x-y) > tol) error stop 'test_nls: '//label
    end subroutine assert_close
end program test_nls
