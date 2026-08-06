program test_durbin
    use ltsa, only : dp, dl_ar_result, ltsa_error, dl_acf_to_ar, dl_residuals, dl_loglikelihood, dl_simulate
    implicit none
    type(dl_ar_result) :: fit
    type(ltsa_error) :: error
    real(dp), allocatable :: r(:), z(:), residuals(:), variances(:), simulated(:)
    real(dp) :: ll, ll_direct, s, phi, gamma0
    real(dp) :: innovations(6)
    integer :: i

    phi = 0.8_dp
    allocate(r(4))
    r = [(phi**i, i=1,4)]
    fit = dl_acf_to_ar(r)
    call assert_true(fit%error%ok(), 'DLAcfToAR failed')
    call assert_close(fit%phi(1), phi, 1.0e-12_dp, 'AR(1) coefficient')
    call assert_max(abs(fit%phi(2:)), 1.0e-12_dp, 'higher AR coefficients')
    call assert_close(fit%pacf(1), phi, 1.0e-12_dp, 'PACF lag one')
    call assert_max(abs(fit%pacf(2:)), 1.0e-12_dp, 'higher PACF')
    call assert_max(abs(fit%prediction_variance-(1.0_dp-phi*phi)), 1.0e-12_dp, 'prediction variances')

    gamma0 = 1.0_dp/(1.0_dp-phi*phi)
    deallocate(r)
    allocate(r(6), z(6))
    r = [(gamma0*phi**(i-1), i=1,6)]
    innovations = [0.25_dp,-0.4_dp,1.2_dp,0.0_dp,-0.7_dp,0.3_dp]
    z(1) = sqrt(gamma0)*innovations(1)
    do i = 2, 6
        z(i) = phi*z(i-1)+innovations(i)
    end do
    call dl_residuals(r,z,residuals,error,standardized=.false.,prediction_variances=variances)
    call assert_true(error%ok(), 'DLResiduals failed')
    call assert_close(residuals(1),z(1),1.0e-11_dp,'initial raw residual')
    call assert_max(abs(residuals(2:)-innovations(2:)), 1.0e-11_dp, 'raw residuals')
    call assert_close(variances(1), gamma0, 1.0e-12_dp, 'initial variance')
    call assert_max(abs(variances(2:)-1.0_dp), 1.0e-12_dp, 'innovation variances')

    ll = dl_loglikelihood(r,z,error)
    call assert_true(error%ok(), 'DLLoglikelihood failed')
    s = innovations(1)**2 + sum(innovations(2:)**2)
    ll_direct = -3.0_dp*log(s/6.0_dp)-0.5_dp*log(gamma0)
    call assert_close(ll,ll_direct,1.0e-11_dp,'concentrated loglikelihood')

    call dl_simulate(6,r,simulated,error,innovations)
    call assert_true(error%ok(), 'DLSimulate failed')
    call assert_max(abs(simulated-z),1.0e-11_dp,'DLSimulate recursion')

    print '(a)', 'test_durbin: PASS'
contains
    subroutine assert_true(condition,message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        if (.not. condition) error stop message
    end subroutine assert_true
    subroutine assert_close(actual,expected,tolerance,message)
        real(dp), intent(in) :: actual,expected,tolerance
        character(len=*), intent(in) :: message
        if (abs(actual-expected) > tolerance) then
            print *, message, actual, expected
            error stop
        end if
    end subroutine assert_close
    subroutine assert_max(values,tolerance,message)
        real(dp), intent(in) :: values(:),tolerance
        character(len=*), intent(in) :: message
        if (size(values) > 0 .and. maxval(abs(values)) > tolerance) then
            print *, message, maxval(abs(values))
            error stop
        end if
    end subroutine assert_max
end program test_durbin
