program test_risk
    use actuar
    implicit none
    type(aggregate_dist_t) :: agg
    type(bstraub_result_t) :: cr
    type(bayes_linear_result_t) :: br
    real(dp) :: fx(2), ratios(3,4), weights(3,4), x(4)
    integer :: i

    fx = [0.0_dp, 1.0_dp]
    agg = panjer_poisson(fx, 2.5_dp, tol=1.0e-12_dp, maxit=80)
    do i=0,min(12,size(agg%pmf)-1)
        call check_close(agg%pmf(i+1), poisson_pmf(i,2.5_dp), 2.0e-11_dp, 'Panjer Poisson')
    end do
    call check_close(agg%mean(),2.5_dp,1.0e-8_dp,'aggregate mean')
    call check_close(agg%variance(),2.5_dp,2.0e-7_dp,'aggregate variance')
    call check_true(aggregate_var(agg,0.95_dp)>=4.0_dp,'aggregate VaR')
    call check_true(aggregate_cte(agg,0.95_dp)>=aggregate_var(agg,0.95_dp),'aggregate CTE')

    ratios = reshape([1.0_dp,1.1_dp,0.9_dp,1.0_dp, 1.3_dp,1.2_dp,1.4_dp,1.3_dp, &
                      0.8_dp,0.9_dp,0.7_dp,0.8_dp], shape(ratios))
    weights = 1.0_dp
    cr = bstraub_fit(ratios,weights,iterative=.true.)
    call check_true(all(cr%credibility>=0.0_dp .and. cr%credibility<=1.0_dp),'credibility range')
    call check_true(size(cr%premium)==3,'credibility premium size')

    x=[2.0_dp,1.0_dp,3.0_dp,2.0_dp]
    br = bayes_poisson_gamma(x,2.0_dp,0.5_dp)
    call check_true(br%premium>0.0_dp,'Bayes premium')

    print '(a)', 'test_risk: PASS'
contains
    subroutine check_close(got, expected, eps, name)
        real(dp),intent(in)::got,expected,eps
        character(*),intent(in)::name
        if(abs(got-expected)>eps*max(1.0_dp,abs(expected))) then
            print *, 'FAIL ',trim(name),got,expected
            error stop 1
        end if
    end subroutine check_close
    subroutine check_true(ok,name)
        logical,intent(in)::ok
        character(*),intent(in)::name
        if(.not.ok) then
            print *, 'FAIL ',trim(name)
            error stop 1
        end if
    end subroutine check_true
end program test_risk
