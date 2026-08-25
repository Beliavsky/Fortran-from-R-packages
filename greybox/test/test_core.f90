program test_core
    use greybox_kinds, only: dp
    use greybox_special, only: normal_pdf
    use greybox_distributions
    use greybox_metrics
    use greybox_utils
    implicit none
    integer :: failures, i
    real(dp) :: p, q, ref
    real(dp) :: a(4), f(4), b(4), low(4), up(4)
    real(dp) :: poly1(2), poly2(3), prod(4)
    real(dp) :: phi(1), beta(2), dm(4)

    failures = 0
    call check_close(dlaplace(0.0_dp,0.0_dp,2.0_dp),0.25_dp,1.0e-13_dp,'laplace density')
    call check_close(plaplace(0.0_dp,0.0_dp,2.0_dp),0.5_dp,1.0e-13_dp,'laplace cdf')
    call check_close(qlaplace(0.75_dp,0.0_dp,2.0_dp),2.0_dp*log(2.0_dp),1.0e-12_dp,'laplace quantile')

    do i=1,5
        p = 0.1_dp + 0.2_dp*real(i-1,dp)
        q = qalaplace(p,1.0_dp,2.0_dp,0.3_dp)
        call check_close(palaplace(q,1.0_dp,2.0_dp,0.3_dp),p,2.0e-11_dp,'alaplace inversion')
    end do

    ref = normal_pdf(0.4_dp,-0.2_dp,1.7_dp/sqrt(2.0_dp))
    call check_close(dgnorm(0.4_dp,-0.2_dp,1.7_dp,2.0_dp),ref,2.0e-13_dp,'gnorm normal special case')
    do i=1,7
        p = 0.05_dp + 0.15_dp*real(i-1,dp)
        q = qgnorm(p,0.2_dp,1.1_dp,1.4_dp)
        call check_close(pgnorm(q,0.2_dp,1.1_dp,1.4_dp),p,3.0e-9_dp,'gnorm inversion')
        q = qs(p,0.1_dp,0.8_dp)
        call check_close(ps(q,0.1_dp,0.8_dp),p,3.0e-9_dp,'S inversion')
    end do

    do i=1,5
        p = 0.1_dp + 0.2_dp*real(i-1,dp)
        q = qfnorm(p,0.7_dp,1.2_dp)
        call check_close(pfnorm(q,0.7_dp,1.2_dp),p,2.0e-10_dp,'folded normal inversion')
        q = qbcnorm(p,0.2_dp,0.7_dp,0.0_dp)
        call check_close(pbcnorm(q,0.2_dp,0.7_dp,0.0_dp),p,2.0e-10_dp,'Box-Cox normal inversion')
        q = qlogitnorm(p,-0.2_dp,0.8_dp)
        call check_close(plogitnorm(q,-0.2_dp,0.8_dp),p,2.0e-10_dp,'logit-normal inversion')
        q = qtplnorm(p,0.1_dp,0.5_dp,3.0_dp)
        call check_close(ptplnorm(q,0.1_dp,0.5_dp,3.0_dp),p,2.0e-10_dp,'three-parameter lognormal inversion')
    end do
    p = prectnorm(0.0_dp,-0.5_dp,1.2_dp)
    call check_close(qrectnorm(0.5_dp*p,-0.5_dp,1.2_dp),0.0_dp,1.0e-14_dp,'rectified point mass')

    a = [1.0_dp,2.0_dp,4.0_dp,8.0_dp]
    f = [1.5_dp,1.0_dp,5.0_dp,7.0_dp]
    b = [1.0_dp,2.5_dp,3.0_dp,9.0_dp]
    low = f-1.0_dp
    up = f+1.0_dp
    call check_close(me(a,f),0.125_dp,1.0e-14_dp,'ME')
    call check_close(mae(a,f),0.875_dp,1.0e-14_dp,'MAE')
    call check_close(mse(a,f),0.8125_dp,1.0e-14_dp,'MSE')
    call check_true(mis(a,low,up,0.1_dp)>0.0_dp,'MIS positive')
    call check_true(gmrae(a,f,b)>0.0_dp,'GMRAE positive')
    call check_true(ham(a,0.0_dp)>0.0_dp,'half absolute moment')

    poly1 = [1.0_dp,2.0_dp]
    poly2 = [3.0_dp,4.0_dp,5.0_dp]
    prod = polyprod(poly1,poly2)
    call check_vec(prod,[3.0_dp,10.0_dp,13.0_dp,10.0_dp],1.0e-14_dp,'polyprod')
    phi = [0.5_dp]
    beta = [1.0_dp,2.0_dp]
    dm = dyn_mult_calc(phi,beta,4)
    call check_vec(dm,[1.0_dp,2.5_dp,1.25_dp,0.625_dp],1.0e-14_dp,'dynamic multipliers')

    if (failures /= 0) then
        write(*,'(a,i0)') 'test_core: FAIL ',failures
        error stop 1
    end if
    write(*,'(a)') 'test_core: PASS'

contains
    subroutine check_close(x,y,tol,name)
        real(dp), intent(in) :: x,y,tol
        character(len=*), intent(in) :: name
        if (abs(x-y) > tol*max(1.0_dp,abs(y))) then
            failures = failures+1
            write(*,'(a,2es20.10)') trim(name)//': ',x,y
        end if
    end subroutine check_close
    subroutine check_vec(x,y,tol,name)
        real(dp), intent(in) :: x(:),y(:),tol
        character(len=*), intent(in) :: name
        if (maxval(abs(x-y)) > tol*max(1.0_dp,maxval(abs(y)))) then
            failures=failures+1
            write(*,'(a)') trim(name)//': mismatch'
        end if
    end subroutine check_vec
    subroutine check_true(ok,name)
        logical,intent(in)::ok
        character(len=*),intent(in)::name
        if(.not.ok)then;failures=failures+1;write(*,'(a)')trim(name)//': failed';end if
    end subroutine check_true
end program test_core
