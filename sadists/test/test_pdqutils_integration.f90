program test_pdqutils_integration
    use sadists
    use pdqutils, only : dapx_edgeworth, papx_edgeworth, qapx_cf, moment2cumulant, cumulant2moment
    implicit none
    real(dp), parameter :: tol = 2.0e-13_dp
    real(dp) :: kappa(6), x, p, q1, q2
    real(dp) :: w(2), df(2), ncp(2), pow(2)
    real(dp) :: moms(5), c1(5), c2(5), m1(5), m2(5)
    integer :: failures

    failures = 0
    w = [0.75_dp, 1.25_dp]
    df = [3.5_dp, 7.0_dp]
    ncp = [0.4_dp, 1.1_dp]
    pow = [1.0_dp, 1.0_dp]
    call sumchisqpow_cumulants(w, df, ncp, pow, kappa)

    x = 9.25_dp
    call check_close(dsumchisqpow(x,w,df,ncp,pow), &
        dapx_edgeworth(x,kappa,support_lo=0.0_dp), tol, 'density', failures)
    call check_close(psumchisqpow(x,w,df,ncp,pow), &
        papx_edgeworth(x,kappa,support_lo=0.0_dp), tol, 'cdf', failures)

    p = 0.73_dp
    q1 = qsumchisqpow(p,w,df,ncp,pow)
    q2 = qapx_cf(p,kappa,support_lo=0.0_dp)
    call check_close(q1,q2,tol,'quantile',failures)

    moms = [1.1_dp, 2.3_dp, 5.7_dp, 16.2_dp, 51.4_dp]
    call moments_to_cumulants(moms,c1)
    c2 = moment2cumulant(moms)
    call check_vec(c1,c2,tol,'moment2cumulant',failures)
    call cumulants_to_moments(c1,m1)
    m2 = cumulant2moment(c1)
    call check_vec(m1,m2,tol,'cumulant2moment',failures)

    if (failures /= 0) then
        print *, 'test_pdqutils_integration: FAIL', failures
        error stop 1
    end if
    print *, 'test_pdqutils_integration: PASS'

contains

    subroutine check_close(a,b,eps,label,failures)
        real(dp), intent(in) :: a,b,eps
        character(*), intent(in) :: label
        integer, intent(inout) :: failures
        if (abs(a-b) > eps*max(1.0_dp,abs(a),abs(b))) then
            print *, 'FAIL ', trim(label), a, b
            failures = failures + 1
        end if
    end subroutine check_close

    subroutine check_vec(a,b,eps,label,failures)
        real(dp), intent(in) :: a(:),b(:),eps
        character(*), intent(in) :: label
        integer, intent(inout) :: failures
        integer :: i
        do i=1,size(a)
            call check_close(a(i),b(i),eps,label,failures)
        end do
    end subroutine check_vec

end program test_pdqutils_integration
